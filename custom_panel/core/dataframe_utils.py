"""
DataFrame utility functions to eliminate DRY violations.

This module provides common DataFrame operations used throughout the codebase
to reduce code duplication and improve maintainability.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import pandas as pd
else:
    import pandas as pd


def safe_column_sum(df: pd.DataFrame, column: str) -> int:
    """
    Safely sum non-NaN values in a column if it exists.

    Args:
        df: DataFrame to check
        column: Column name to sum

    Returns:
        Sum of non-NaN values, or 0 if column doesn't exist
    """
    return (~df[column].isna()).sum() if column in df.columns else 0


def safe_column_count(df: pd.DataFrame, column: str) -> int:
    """
    Safely count non-NaN values in a column if it exists.

    Args:
        df: DataFrame to check
        column: Column name to count

    Returns:
        Count of non-NaN values, or 0 if column doesn't exist
    """
    return df[column].notna().sum() if column in df.columns else 0


def safe_bool_count(df: pd.DataFrame, column: str) -> int:
    """
    Safely count True values in a boolean column if it exists.

    Args:
        df: DataFrame to check
        column: Column name to count

    Returns:
        Count of True values, or 0 if column doesn't exist
    """
    if column in df.columns:
        return (df[column] == True).sum()  # noqa: E712
    return 0


def safe_column_value(df: pd.DataFrame, column: str, default: Any = None) -> Any:
    """
    Safely get column values with default if column doesn't exist.

    Args:
        df: DataFrame to check
        column: Column name to get
        default: Default value if column doesn't exist

    Returns:
        Column values or default
    """
    return df[column] if column in df.columns else default


def safe_column_mean(df: pd.DataFrame, column: str, default: float = 0.0) -> float:
    """
    Safely calculate mean of a column if it exists.

    Args:
        df: DataFrame to check
        column: Column name to calculate mean for
        default: Default value if column doesn't exist or has no valid values

    Returns:
        Mean value or default
    """
    if column in df.columns and not df[column].empty:
        mean_val = df[column].mean()
        return mean_val if not pd.isna(mean_val) else default
    return default


def safe_column_max(df: pd.DataFrame, column: str, default: Any = None) -> Any:
    """
    Safely get maximum value of a column if it exists.

    Args:
        df: DataFrame to check
        column: Column name to get max for
        default: Default value if column doesn't exist or has no valid values

    Returns:
        Maximum value or default
    """
    if column in df.columns and not df[column].empty:
        max_val = df[column].max()
        return max_val if not pd.isna(max_val) else default
    return default


def safe_nan_check(
    value: Any, converter_func: Callable[[Any], Any] | None = None, default: Any = None
) -> Any:
    """
    Safely handle NaN values with optional conversion.

    Args:
        value: Value to check
        converter_func: Optional function to convert non-NaN values
        default: Default value for NaN values

    Returns:
        Converted value, original value, or default
    """
    if pd.isna(value):
        return default
    return converter_func(value) if converter_func else value


def safe_int_conversion(value: Any, default: int | None = None) -> int | None:
    """
    Safely convert value to int handling NaN.

    Args:
        value: Value to convert
        default: Default value for NaN values

    Returns:
        Integer value or default
    """
    return safe_nan_check(value, int, default)


def safe_float_conversion(value: Any, default: float | None = None) -> float | None:
    """
    Safely convert value to float handling NaN.

    Args:
        value: Value to convert
        default: Default value for NaN values

    Returns:
        Float value or default
    """
    return safe_nan_check(value, float, default)


def safe_str_conversion(value: Any, default: str = "") -> str:
    """
    Safely convert value to string handling NaN.

    Args:
        value: Value to convert
        default: Default value for NaN values

    Returns:
        String value or default
    """
    return safe_nan_check(value, str, default)


def process_dataframe_rows(
    df: pd.DataFrame,
    processor_func: Callable[[pd.Series[Any]], Any],
    filter_func: Callable[[pd.Series[Any]], bool] | None = None,
) -> list[Any]:
    """
    Abstract DataFrame row iteration with optional filtering.

    Args:
        df: DataFrame to process
        processor_func: Function to apply to each row
        filter_func: Optional function to filter rows

    Returns:
        List of processed results
    """
    results = []
    for _, row in df.iterrows():
        if filter_func is None or filter_func(row):
            result = processor_func(row)
            if result is not None:
                results.append(result)
    return results


def extract_numeric_data(
    df: pd.DataFrame, columns: list[str]
) -> dict[str, list[float]]:
    """
    Extract numeric data from specified columns, filtering out NaN values.

    Args:
        df: DataFrame to extract from
        columns: List of column names to extract

    Returns:
        Dictionary mapping column names to lists of numeric values
    """
    result = {}
    for column in columns:
        if column in df.columns:
            numeric_values = [
                float(value) for value in df[column] if not pd.isna(value)
            ]
            result[column] = numeric_values
        else:
            result[column] = []
    return result


def get_column_statistics(df: pd.DataFrame, column: str) -> dict[str, Any]:
    """
    Get comprehensive statistics for a column.

    Args:
        df: DataFrame to analyze
        column: Column name to get statistics for

    Returns:
        Dictionary with count, mean, max, min statistics
    """
    if column not in df.columns:
        return {"count": 0, "mean": None, "max": None, "min": None}

    non_na_values = df[column].dropna()
    if non_na_values.empty:
        return {"count": 0, "mean": None, "max": None, "min": None}

    return {
        "count": len(non_na_values),
        "mean": non_na_values.mean() if non_na_values.dtype.kind in "biufc" else None,
        "max": non_na_values.max(),
        "min": non_na_values.min(),
    }


def check_columns_exist(
    df: pd.DataFrame, required_columns: list[str]
) -> dict[str, bool]:
    """
    Check which columns exist in the DataFrame.

    Args:
        df: DataFrame to check
        required_columns: List of column names to check for

    Returns:
        Dictionary mapping column names to existence boolean
    """
    return {col: col in df.columns for col in required_columns}


def filter_existing_columns(df: pd.DataFrame, columns: list[str]) -> list[str]:
    """
    Filter list of columns to only those that exist in the DataFrame.

    Args:
        df: DataFrame to check
        columns: List of column names to filter

    Returns:
        List of column names that exist in the DataFrame
    """
    return [col for col in columns if col in df.columns]
