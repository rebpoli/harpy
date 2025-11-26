
def format_label(t_seconds):
    """
    Format time in seconds as a nice label with appropriate units.

    Parameters:
    -----------
    t_seconds : float
        Time in seconds

    Returns:
    --------
    str : Formatted label like "5.2 hours", "3.1 days", "2.5 months", "1.2 years"
    """
    # Time conversion constants
    HOUR = 3600
    DAY = 86400
    MONTH = 2592000  # 30 days
    YEAR = 31536000  # 365 days

    # Choose the most appropriate unit
    if t_seconds < 2 * DAY:  # Less than 2 days → use hours
        value = t_seconds / HOUR
        if value < 1:
            unit = "hour"
            precision = 2
        else:
            unit = "hours"
            precision = 1

    elif t_seconds < 60 * DAY:  # Less than 60 days → use days
        value = t_seconds / DAY
        unit = "days" if value >= 2 else "day"
        precision = 1

    elif t_seconds < 2 * YEAR:  # Less than 2 years → use months
        value = t_seconds / MONTH
        unit = "months" if value >= 2 else "month"
        precision = 1

    else:  # 2 years or more → use years
        value = t_seconds / YEAR
        unit = "years" if value >= 2 else "year"
        precision = 1 if value < 10 else 0

    # Format the value
    if precision == 0:
        return f"t = {value:.0f} {unit}"
    elif precision == 1:
        return f"t = {value:.1f} {unit}"
    else:
        return f"t = {value:.2f} {unit}"
