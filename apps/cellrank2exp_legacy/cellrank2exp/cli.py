from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def _find_project_root(start: Path) -> Path:
    start = start.resolve()
    for candidate in (start, *start.parents):
        if (candidate / "package.json").is_file():
            return candidate
    raise FileNotFoundError(
        "Could not find project root (missing package.json). "
        "Run this from the repo root."
    )


def _require_cmd(cmd: str, install_hint: str) -> None:
    if shutil.which(cmd) is None:
        raise SystemExit(f"Missing `{cmd}` on PATH. {install_hint}")


def _run(cmd: list[str], *, cwd: Path) -> None:
    subprocess.run(cmd, cwd=str(cwd), check=True)


def _maybe_warn_missing_api_key(project_root: Path) -> None:
    env_file = project_root / ".env.local"
    if not env_file.is_file():
        return
    try:
        content = env_file.read_text(encoding="utf-8")
    except OSError:
        return
    if "GEMINI_API_KEY=PLACEHOLDER_API_KEY" in content:
        print(
            "Warning: `.env.local` still contains `GEMINI_API_KEY=PLACEHOLDER_API_KEY`.",
            file=sys.stderr,
        )


def _ensure_node_modules(project_root: Path) -> None:
    if (project_root / "node_modules").is_dir():
        return
    print("Installing npm dependencies…", file=sys.stderr)
    _run(["npm", "install"], cwd=project_root)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="cellrank2exp")
    sub = parser.add_subparsers(dest="command", required=True)

    install_parser = sub.add_parser("install", help="Install npm dependencies")
    install_parser.add_argument(
        "--cwd",
        default=None,
        help="Project directory (defaults to current directory / parents).",
    )

    start_parser = sub.add_parser("start", help="Start the Vite dev server")
    start_parser.add_argument(
        "--cwd",
        default=None,
        help="Project directory (defaults to current directory / parents).",
    )
    start_parser.add_argument("--host", default=None, help="Vite host (passed through).")
    start_parser.add_argument("--port", type=int, default=None, help="Vite port.")
    start_parser.add_argument(
        "--no-install",
        action="store_true",
        help="Do not auto-run `npm install` when node_modules is missing.",
    )
    start_parser.add_argument(
        "vite_args",
        nargs=argparse.REMAINDER,
        help="Extra args for Vite (use `--` before these args).",
    )

    args = parser.parse_args(argv)

    start_dir = Path(args.cwd).expanduser() if getattr(args, "cwd", None) else Path.cwd()
    project_root = _find_project_root(start_dir)

    _require_cmd(
        "npm",
        "Install Node.js/npm or recreate the conda env from `environment.yml`.",
    )

    if args.command == "install":
        _ensure_node_modules(project_root)
        return

    if args.command == "start":
        _maybe_warn_missing_api_key(project_root)
        if not args.no_install:
            _ensure_node_modules(project_root)

        cmd = ["npm", "run", "dev", "--"]
        if args.host:
            cmd += ["--host", args.host]
        if args.port:
            cmd += ["--port", str(args.port)]
        passthrough = list(args.vite_args or [])
        if passthrough[:1] == ["--"]:
            passthrough = passthrough[1:]
        cmd += passthrough

        _run(cmd, cwd=project_root)
        return

    raise SystemExit(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
