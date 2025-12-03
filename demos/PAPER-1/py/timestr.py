
#  
#
#
def format_time_duration(seconds):
    if seconds < 0:
        return f"-{format_time_duration(-seconds)}"

    # Time units in seconds
    YEAR = 365.25 * 24 * 3600    # ~31,557,600 seconds (accounting for leap years)
    MONTH = 30.44 * 24 * 3600    # ~2,629,440 seconds (average month)
    DAY = 24 * 3600              # 86,400 seconds
    HOUR = 3600                  # 3,600 seconds
    MINUTE = 60                  # 60 seconds

    # Convert to different units
    years = int(seconds // YEAR)
    seconds %= YEAR

    months = int(seconds // MONTH)
    seconds %= MONTH

    days = int(seconds // DAY)
    seconds %= DAY

    hours = int(seconds // HOUR)
    seconds %= HOUR

    minutes = int(seconds // MINUTE)
    remaining_seconds = int(seconds % MINUTE)

    # Build formatted string
    parts = []

    if years > 0:
        parts.append(f"{years}y")
    if months > 0:
        parts.append(f"{months}m")
    if days > 0:
        parts.append(f"{days}d")
    if hours > 0:
        parts.append(f"{hours}h")
    if minutes > 0:
        parts.append(f"{minutes}min")
    if remaining_seconds > 0 or len(parts) == 0:  # Always show seconds if nothing else
        parts.append(f"{remaining_seconds}s")

    # Return appropriate level of detail
    if years > 0:
        return " ".join(parts[:3])  # Show years, months, days
    elif months > 0:
        return " ".join(parts[:3])  # Show months, days, hours
    elif days > 0:
        return " ".join(parts[:3])  # Show days, hours, minutes
    elif hours > 0:
        return " ".join(parts[:2])  # Show hours, minutes
    else:
        return " ".join(parts[:2])  # Show minutes, seconds

