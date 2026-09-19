import { Elysia } from "elysia";

const dir = import.meta.dir;
const root = `${dir}/..`;
const python = `${root}/.venv/bin/python`;
const port = Number(Bun.env.PORT ?? 3100);

// Fixed whitelist. The route never builds a command from request input, so a
// request cannot run anything but these.
const EXPERIMENTS: Record<string, string> = {
  e1: "experiments/e1_regimes.py",
  e2: "experiments/e2_threshold.py",
  e3: "experiments/e3_snr_vs_noise.py",
  e4: "experiments/e4_robustness.py",
  e5: "experiments/e5_parameter_tuned.py",
  data: "web/export_data.py",
};

type Job = { name: string; log: string; started: number; done: boolean; code: number | null };
let job: Job | null = null;

const file = (name: string, type: string) =>
  new Response(Bun.file(`${dir}/${name}`), {
    headers: { "content-type": type, "cache-control": "no-store" },
  });

async function drain(job: Job, stream: ReadableStream<Uint8Array> | null) {
  if (!stream) return;
  const decoder = new TextDecoder();
  for await (const chunk of stream) {
    job.log = (job.log + decoder.decode(chunk)).slice(-4000);
  }
}

async function runScript(job: Job, script: string) {
  const proc = Bun.spawn([python, "-u", script], {
    cwd: root,
    stdout: "pipe",
    stderr: "pipe",
  });
  await Promise.all([drain(job, proc.stdout), drain(job, proc.stderr)]);
  return await proc.exited;
}

async function start(name: string) {
  const current: Job = { name, log: "", started: Date.now(), done: false, code: null };
  job = current;
  let code = await runScript(current, EXPERIMENTS[name]);
  // Refresh the dashboard numbers after any experiment.
  if (code === 0 && name !== "data") {
    current.log += "\n";
    code = await runScript(current, EXPERIMENTS.data);
  }
  current.code = code;
  current.done = true;
}

new Elysia()
  .get("/", () => file("index.html", "text/html; charset=utf-8"))
  .get("/data.js", () => file("data.js", "text/javascript; charset=utf-8"))
  .get("/figures/:name", ({ params, set }) => {
    if (!/^[a-z0-9_]+\.png$/.test(params.name)) {
      set.status = 404;
      return "not found";
    }
    return new Response(Bun.file(`${root}/figures/${params.name}`), {
      headers: { "content-type": "image/png", "cache-control": "no-store" },
    });
  })
  .get("/run/status", () => ({
    running: job ? !job.done : false,
    job: job && {
      name: job.name,
      seconds: Math.round((Date.now() - job.started) / 1000),
      done: job.done,
      code: job.code,
      log: job.log.split("\n").slice(-12).join("\n"),
    },
    experiments: Object.keys(EXPERIMENTS),
  }))
  .post("/run/:name", ({ params, set }) => {
    if (!(params.name in EXPERIMENTS)) {
      set.status = 404;
      return { error: `unknown target ${params.name}` };
    }
    if (job && !job.done) {
      set.status = 409;
      return { error: `${job.name} is still running` };
    }
    start(params.name);
    return { started: params.name };
  })
  .listen(port);

console.log(`dashboard on http://localhost:${port}`);
