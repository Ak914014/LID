from pydantic_settings import BaseSettings

class Settings(BaseSettings):
    APP_NAME: str = "Pocket Diagram PoC"
    ENV: str = "dev"
    ALLOWED_ORIGINS: list[str] = [
        "http://localhost:5173",
        "http://127.0.0.1:5173",
    ]
    UPLOAD_DIR: str = "storage/uploads"

settings = Settings()