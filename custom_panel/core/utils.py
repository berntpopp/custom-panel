"""
Utility functions for the custom-panel tool.

This module provides common utility functions used across the application.
"""

import math
from typing import Any


def normalize_count(
    count: int, method: str = "logistic", params: dict[str, Any] | None = None
) -> float:
    """
    Normalize a count value to a confidence score between 0 and 1.

    Args:
        count: The count value to normalize (e.g., number of sources)
        method: The normalization method ("logistic" or "linear")
        params: Parameters for the normalization method
            - For "logistic": requires "k" (steepness) and "x0" (midpoint)
            - For "linear": requires "max_count" (count for full confidence)

    Returns:
        Normalized confidence score between 0 and 1

    Raises:
        ValueError: If method is unknown or required parameters are missing
    """
    if params is None:
        params = {}

    if count < 0:
        raise ValueError(f"Count must be non-negative, got {count}")

    if method == "logistic":
        # Logistic normalization: 1 / (1 + e^(-k*(x-x0)))
        k = params.get("k")
        x0 = params.get("x0")

        if k is None or x0 is None:
            raise ValueError("Logistic normalization requires 'k' and 'x0' parameters")

        try:
            # Prevent overflow in exponential calculation
            exponent = -k * (count - x0)
            if exponent > 100:  # e^100 is effectively infinity
                return 0.0
            elif exponent < -100:  # e^-100 is effectively 0
                return 1.0
            else:
                return 1 / (1 + math.exp(exponent))
        except (OverflowError, ZeroDivisionError):
            # Handle edge cases
            return 1.0 if count > x0 else 0.0

    elif method == "linear":
        # Linear normalization: min(count / max_count, 1.0)
        max_count = params.get("max_count")

        if max_count is None:
            raise ValueError("Linear normalization requires 'max_count' parameter")

        if max_count <= 0:
            raise ValueError(f"max_count must be positive, got {max_count}")

        return min(count / max_count, 1.0)

    else:
        raise ValueError(f"Unknown normalization method: {method}")


def get_source_group_from_config(
    source_name: str, config: dict[str, Any]
) -> str | None:
    """
    Get the source group name for a given source from the configuration.

    Args:
        source_name: The name of the specific source
        config: The configuration dictionary

    Returns:
        The source group name if found, None otherwise
    """
    data_sources = config.get("data_sources", {})

    # Check each source group
    for group_name, group_config in data_sources.items():
        if not isinstance(group_config, dict):
            continue

        # Check if this is a source group
        if group_config.get("source_group", False):
            # Check if the source_name is in the panels list
            panels = group_config.get("panels", [])
            for panel in panels:
                if isinstance(panel, dict) and panel.get("name") == source_name:
                    return group_name
        else:
            # For standalone sources, the group name might match the source name
            # or it might have a specific name field
            if group_name == source_name or group_config.get("name") == source_name:
                return group_name

    return None
