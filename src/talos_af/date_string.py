"""
This is a placeholder, completely base class to prevent circular imports
"""

import zoneinfo
from datetime import datetime
from functools import cache


@cache
def get_date_string():
    """Cached date string getter."""
    return datetime.now(tz=zoneinfo.ZoneInfo('Australia/Brisbane')).strftime('%Y-%m-%d')
