import sys
import logging
from typing import Optional

_logger: Optional[logging.Logger] = None

def setup_logger(
        name: str = "synnet",
        level: str = "INFO",
        log_file: Optional[str] = None,
        fmt: str = "%(asctime)s [%(levelname)s] %(message)s",
        datefmt: str = "%H:%M:%S",
) -> logging.Logger:
    global _logger
    if _logger is not None:
        return _logger

    # 创建日志器
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))

    # 避免重复添加处理器
    if logger.handlers:
        return logger

    # 控制台处理器
    console = logging.StreamHandler(sys.stderr)
    console.setLevel(logging.DEBUG)
    console.setFormatter(logging.Formatter(fmt, datefmt))
    logger.addHandler(console)

    # 文件处理器 (可选)
    if log_file:
        from pathlib import Path
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, mode='a', encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(fmt.replace('%(asctime)s', '%(filename)s:%(lineno)d')))
        logger.addHandler(file_handler)

    _logger = logger
    return logger

def get_logger(name: str = None) -> logging.Logger:
    return setup_logger(name=name) if name else (_logger or setup_logger())

def debug(msg, *args, **kwargs):
    get_logger().debug(msg, *args, **kwargs)

def info(msg, *args, **kwargs):
    get_logger().info(msg, *args, **kwargs)

def warning(msg, *args, **kwargs):
    get_logger().warning(msg, *args, **kwargs)

def error(msg, *args, **kwargs):
    get_logger().error(msg, *args, **kwargs)

def success(msg, *args, **kwargs):
    info(f"{msg}", *args, **kwargs)

if __name__ == "__main__":
    logger = setup_logger(level="DEBUG")
    debug("This is DEBUG")
    info("This is INFO")
    warning("This is WARNING")
    error("This is ERROR")
    success("This is SUCCESS")