import time
from fastapi import Request
from starlette.middleware.base import BaseHTTPMiddleware

class RequestTimingMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request: Request, call_next):
        start = time.time()
        resp = await call_next(request)
        resp.headers["X-Process-Time-ms"] = f"{(time.time() - start)*1000:.2f}"
        return resp