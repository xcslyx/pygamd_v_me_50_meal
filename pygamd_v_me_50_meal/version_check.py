# pygamd_v_me_50_meal/version_check.py

import importlib.metadata
import json
import os
import sys
import urllib.request
from datetime import datetime, timezone
from packaging.version import Version

OWNER = "lyxlyxlyxxx"
REPO = "pygamd_v_me_50_meal"
PKG_NAME = "pygamd_v_me_50_meal"

def _debug(msg: str):
    if os.getenv("PYGAMD_UPDATE_DEBUG"):
        print(f"[pygamd:update] {msg}", file=sys.stderr)


def _get_local_version():
    try:
        v = Version(importlib.metadata.version(PKG_NAME))
        _debug(f"Local version detected: {v}")
        return v
    except Exception as e:
        _debug(f"Failed to read local version: {e}")
        return None


def _get_repo_info(timeout):
    url = f"https://gitee.com/api/v5/repos/{OWNER}/{REPO}"
    _debug(f"Querying repo info: {url}")
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        return json.load(resp)


def _get_latest_release_version(timeout):
    url = (
        f"https://gitee.com/api/v5/repos/"
        f"{OWNER}/{REPO}/releases/latest"
    )
    _debug("Querying latest release info")
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        data = json.load(resp)

    tag = data.get("tag_name")
    if not tag:
        _debug("No tag_name found in latest release")
        return None

    try:
        v = Version(tag.lstrip("v"))
        _debug(f"Latest release version: {v}")
        return v
    except Exception as e:
        _debug(f"Invalid tag format '{tag}': {e}")
        return None


def _parse_utc(ts: str):
    return datetime.fromisoformat(ts.replace("Z", "+00:00"))


def check_update(timeout: float = 2.0):
    if os.getenv("PYGAMD_NO_UPDATE_CHECK"):
        _debug("Update check disabled by environment variable")
        return

    try:
        local_version = _get_local_version()
        repo_info = _get_repo_info(timeout)

        # ===== 优先级 1：release / tag =====
        try:
            latest_release = _get_latest_release_version(timeout)
        except Exception as e:
            latest_release = None
            _debug(f"Release query failed: {e}")

        if local_version and latest_release:
            if latest_release > local_version:
                print(
                    f"[pygamd] New version available: "
                    f"{local_version} → {latest_release}\n"
                    f"https://gitee.com/{OWNER}/{REPO}"
                )
            else:
                _debug("Local version is up to date with latest release")
            return

        _debug("Release-based check not available, using pushed_at fallback")

        # ===== fallback：pushed_at =====
        pushed_at = repo_info.get("pushed_at")
        if not pushed_at:
            _debug("No pushed_at field in repo info")
            return

        remote_time = _parse_utc(pushed_at)
        _debug(f"Repository last pushed at: {remote_time.isoformat()}")

        try:
            dist = importlib.metadata.distribution(PKG_NAME)
            local_time = datetime.fromtimestamp(
                dist._path.stat().st_mtime, tz=timezone.utc
            )
            _debug(f"Local install time: {local_time.isoformat()}")
        except Exception as e:
            _debug(f"Failed to determine local install time: {e}")
            return

        if remote_time > local_time:
            print(
                "[pygamd] Repository has been updated since "
                "your local installation.\n"
                f"Please run 'pip install git+https://gitee.com/{OWNER}/{REPO}'"
            )
        else:
            _debug("Local installation is newer than repo push time")

    except Exception as e:
        _debug(f"Unexpected failure during update check: {e}")


if __name__ == "__main__":
    # os.environ["PYGAMD_UPDATE_DEBUG"] = "1"
    check_update()