"""Configuration helpers for the cleaned repository."""

from pathlib import Path
from typing import Any


def find_repo_root(start: Path | None = None) -> Path:
    """Return the repository root containing configs/config.yaml."""
    current = (start or Path.cwd()).resolve()
    for candidate in (current, *current.parents):
        if (candidate / "configs" / "config.yaml").exists():
            return candidate
    raise FileNotFoundError("Could not find repository root containing configs/config.yaml.")


def load_config(config_path: str | Path | None = None) -> dict[str, Any]:
    """Load configs/config.yaml."""
    try:
        import yaml
    except ImportError as exc:
        raise ImportError("PyYAML is required to read configs/config.yaml.") from exc

    path = Path(config_path) if config_path else find_repo_root() / "configs" / "config.yaml"
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def get_nested(config: dict[str, Any], dotted_key: str, default: Any = None) -> Any:
    """Read nested config values using dotted keys such as paths.result_dir."""
    value: Any = config
    for part in dotted_key.split("."):
        if not isinstance(value, dict) or part not in value:
            return default
        value = value[part]
    return value


def resolve_path(config: dict[str, Any], dotted_key: str, default: str | None = None) -> Path:
    """Resolve a config path relative to the repository root."""
    value = get_nested(config, dotted_key, default)
    if value is None:
        raise KeyError(f"Missing path config: {dotted_key}")
    path = Path(value)
    if path.is_absolute():
        return path
    return find_repo_root() / path


def ensure_dir(path: str | Path) -> Path:
    """Create and return an output directory."""
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out

