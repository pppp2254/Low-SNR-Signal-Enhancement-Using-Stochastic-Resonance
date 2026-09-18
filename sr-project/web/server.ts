import { Elysia } from "elysia";

const dir = import.meta.dir;
const file = (name: string, type: string) =>
  new Response(Bun.file(`${dir}/${name}`), {
    headers: { "content-type": type, "cache-control": "no-store" },
  });

const port = Number(Bun.env.PORT ?? 3100);

new Elysia()
  .get("/", () => file("index.html", "text/html; charset=utf-8"))
  .get("/data.js", () => file("data.js", "text/javascript; charset=utf-8"))
  .listen(port);

console.log(`dashboard on http://localhost:${port}`);
