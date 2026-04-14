import os
from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles

from .core.config import settings
from .core.cors import add_cors
from .core.middleware import RequestTimingMiddleware

from .routes.health import router as health_router
from .routes.diagram import router as diagram_router
from .routes.lid import router as lid_router
from .routes.viz3d import router as viz3d_router

def create_app():
    app = FastAPI(title=settings.APP_NAME)

    # CORS
    add_cors(app)

    # Middleware
    app.add_middleware(RequestTimingMiddleware)

    # Routes
    app.include_router(health_router)
    app.include_router(diagram_router)
    app.include_router(lid_router)
    app.include_router(viz3d_router, prefix="/api")

    # Static mount for uploaded PDBs (NGL will load from here)
    os.makedirs(settings.UPLOAD_DIR, exist_ok=True)
    app.mount("/static", StaticFiles(directory=settings.UPLOAD_DIR), name="static")

    return app

app = create_app()