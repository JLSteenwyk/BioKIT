import json
from typing import Any


class BaseService(object):
    VALID_OUTPUT_FORMATS = ("tsv", "json", "yaml")

    def __init__(self) -> None:
        pass

    def process_args(self, args: object) -> dict:
        raise NotImplementedError()

    def run(self) -> None:
        raise NotImplementedError()

    @classmethod
    def normalize_output_format(cls, output_format: str | None) -> str:
        if output_format is None:
            return "tsv"
        if output_format not in cls.VALID_OUTPUT_FORMATS:
            valid = ", ".join(cls.VALID_OUTPUT_FORMATS)
            raise ValueError(f"Invalid output format {output_format!r}. Expected one of: {valid}")
        return output_format

    @classmethod
    def format_rows(
        cls,
        rows: list[dict[str, Any]],
        output_format: str | None,
        field_order: list[str] | None = None,
    ) -> str:
        output_format = cls.normalize_output_format(output_format)
        if output_format == "json":
            return json.dumps(rows, indent=2)
        if output_format == "yaml":
            return cls._to_yaml(rows)
        if not rows:
            return ""
        fields = field_order if field_order is not None else list(rows[0].keys())
        lines = ["\t".join(fields)]
        for row in rows:
            lines.append("\t".join(str(row.get(field, "")) for field in fields))
        return "\n".join(lines)

    @classmethod
    def format_object(cls, obj: dict[str, Any], output_format: str | None) -> str:
        output_format = cls.normalize_output_format(output_format)
        if output_format == "json":
            return json.dumps(obj, indent=2)
        if output_format == "yaml":
            return cls._to_yaml(obj)
        lines = [f"{value}\t{label}" for label, value in obj.items()]
        return "\n".join(lines)

    @classmethod
    def _to_yaml(cls, value: Any, indent: int = 0) -> str:
        prefix = " " * indent
        if isinstance(value, dict):
            lines: list[str] = []
            for key, val in value.items():
                if isinstance(val, (dict, list)):
                    lines.append(f"{prefix}{key}:")
                    lines.append(cls._to_yaml(val, indent + 2))
                else:
                    lines.append(f"{prefix}{key}: {cls._yaml_scalar(val)}")
            return "\n".join(lines)
        if isinstance(value, list):
            lines = []
            for item in value:
                if isinstance(item, dict):
                    keys = list(item.keys())
                    if not keys:
                        lines.append(f"{prefix}- {{}}")
                        continue
                    first = keys[0]
                    first_val = item[first]
                    if isinstance(first_val, (dict, list)):
                        lines.append(f"{prefix}- {first}:")
                        lines.append(cls._to_yaml(first_val, indent + 4))
                    else:
                        lines.append(f"{prefix}- {first}: {cls._yaml_scalar(first_val)}")
                    for key in keys[1:]:
                        val = item[key]
                        if isinstance(val, (dict, list)):
                            lines.append(f"{prefix}  {key}:")
                            lines.append(cls._to_yaml(val, indent + 4))
                        else:
                            lines.append(f"{prefix}  {key}: {cls._yaml_scalar(val)}")
                elif isinstance(item, list):
                    lines.append(f"{prefix}-")
                    lines.append(cls._to_yaml(item, indent + 2))
                else:
                    lines.append(f"{prefix}- {cls._yaml_scalar(item)}")
            return "\n".join(lines)
        return f"{prefix}{cls._yaml_scalar(value)}"

    @staticmethod
    def _yaml_scalar(value: Any) -> str:
        if isinstance(value, bool):
            return "true" if value else "false"
        if value is None:
            return "null"
        if isinstance(value, str):
            if value == "":
                return '""'
            if any(char in value for char in [":", "#", "-", "{", "}", "[", "]", "\t", "\n"]):
                return json.dumps(value)
            return value
        return str(value)
