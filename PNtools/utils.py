import pymaid
from functools import wraps

def has_remote_instance(function):
    """Decorator to exit function early if not CatmaidInstance."""
    @wraps(function)
    def wrapper(*args, **kwargs):
        # Get remote instance
        pymaid.set_loggers('ERROR')
        rm = pymaid.utils._eval_remote_instance(None, raise_error=False)
        pymaid.set_loggers('INFO')
        if not rm:
            print('No global CatmaidInstance set. Please define and rerun.')
            return

        # Execute function
        res = function(*args, **kwargs)
        return res
    return wrapper

		
