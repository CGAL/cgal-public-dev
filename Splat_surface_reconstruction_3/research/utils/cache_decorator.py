import os
import pickle
import hashlib
from functools import wraps

def reset_call_counts(cache_dir="../cache"):
    """
    Reset the call counts for all decorated functions.
    """
    calls_file = os.path.join(cache_dir, "_call_counts.pkl")
    if os.path.exists(calls_file):
        os.remove(calls_file)

reset_call_counts()



def disk_cache(cache_dir="../cache", version=1, max_calls=5, enabled=True):
    """
    Cache a function's return value to disk based on its arguments.

    If max_calls is specified, keep track of how many times each decorated
    function is called in a separate file. Raise an exception if the function
    is called more than max_calls times.
    """

    os.makedirs(cache_dir, exist_ok=True)

    # Separate file from the actual cache
    calls_file = os.path.join(cache_dir, "_call_counts.pkl")

    def decorator(func):

        @wraps(func)
        def wrapper(*args, **kwargs):

            # ---------------------------------------------------------
            # Call count
            # ---------------------------------------------------------
            if max_calls is not None:
                if os.path.exists(calls_file):
                    with open(calls_file, "rb") as f:
                        call_counts = pickle.load(f)
                else:
                    call_counts = {}

                function_name = f"{func.__module__}.{func.__qualname__}"
                call_counts[function_name] = call_counts.get(function_name, 0) + 1

                with open(calls_file, "wb") as f:
                    pickle.dump(
                        call_counts,
                        f,
                        protocol=pickle.HIGHEST_PROTOCOL
                    )

                if call_counts[function_name] > max_calls:
                    raise RuntimeError(
                        f"Function '{function_name}' has been called "
                        f"{call_counts[function_name]} times "
                        f"(max_calls={max_calls})."
                    )

            # ---------------------------------------------------------
            # Normal cache logic
            # ---------------------------------------------------------
            key_data = (
                func.__name__,
                version,
                args,
                kwargs,
            )

            key = hashlib.sha256(
                pickle.dumps(
                    key_data,
                    protocol=pickle.HIGHEST_PROTOCOL
                )
            ).hexdigest()

            cache_file = os.path.join(cache_dir, f"{key}.pkl")

            if os.path.exists(cache_file) and enabled:
                print(f"Loading cached result for {func.__name__} from {cache_file}")

                with open(cache_file, "rb") as f:
                    return pickle.load(f)

            result = func(*args, **kwargs)

            with open(cache_file, "wb") as f:
                pickle.dump(
                    result,
                    f,
                    protocol=pickle.HIGHEST_PROTOCOL
                )

            print(f"Saved cache for {func.__name__} to {cache_file}")

            return result

        return wrapper

    return decorator