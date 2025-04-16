import time
import logging
from functools import wraps
from typing import Callable, TypeVar, Any

log = logging.getLogger(__name__)

T = TypeVar('T')

def retry(
    max_retries: int = 3,
    initial_delay: float = 1.0,
    max_delay: float = 10.0,
    backoff_factor: float = 2.0,
    exceptions: tuple = (Exception,)
) -> Callable[[Callable[..., T]], Callable[..., T]]:
    """
    Retry decorator with exponential backoff.
    
    Args:
        max_retries: Maximum number of retries
        initial_delay: Initial delay between retries in seconds
        max_delay: Maximum delay between retries in seconds
        backoff_factor: Factor by which the delay increases after each retry
        exceptions: Tuple of exceptions to catch and retry on
    """
    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> T:
            delay = initial_delay
            last_exception = None
            
            for attempt in range(max_retries + 1):
                try:
                    return func(*args, **kwargs)
                except exceptions as e:
                    last_exception = e
                    if attempt == max_retries:
                        log.error(f"Max retries ({max_retries}) exceeded for {func.__name__}: {str(e)}")
                        raise
                    
                    log.warning(f"Attempt {attempt + 1}/{max_retries} failed for {func.__name__}: {str(e)}")
                    time.sleep(min(delay, max_delay))
                    delay *= backoff_factor
            
            raise last_exception  # This should never be reached due to the raise in the loop
        return wrapper
    return decorator 

# Testing
# if __name__ == "__main__":
#     @retry(max_retries=3, initial_delay=1.0, max_delay=10.0, backoff_factor=2.0)
#     def test_function():
#         print("Test function called")
#         raise Exception("Test exception")
#     test_function()