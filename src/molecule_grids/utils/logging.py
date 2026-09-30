"""Package logger: stdlib logging with Rich output and a ``success`` level."""

import logging

from rich.logging import RichHandler

SUCCESS = 25
logging.addLevelName(SUCCESS, "SUCCESS")


class Logger:
    """Singleton wrapper around the ``molecule_grids`` logger."""

    def __init__(self):
        self._logger = logging.getLogger("molecule_grids")
        if not self._logger.handlers:
            handler = RichHandler(show_path=False, markup=False, rich_tracebacks=True)
            handler.setFormatter(logging.Formatter("%(message)s", datefmt="[%X]"))
            self._logger.addHandler(handler)
        self._logger.setLevel(logging.INFO)
        self._logger.propagate = False

    def set_verbosity(self, verbose):
        """Show info messages when ``verbose`` is true, else warnings and errors only."""
        self._logger.setLevel(logging.INFO if verbose else logging.WARNING)

    def debug(self, msg, *args):
        self._logger.debug(msg, *args)

    def info(self, msg, *args):
        self._logger.info(msg, *args)

    def success(self, msg, *args):
        self._logger.log(SUCCESS, msg, *args)

    def warning(self, msg, *args):
        self._logger.warning(msg, *args)

    def error(self, msg, *args):
        self._logger.error(msg, *args)


logger = Logger()
