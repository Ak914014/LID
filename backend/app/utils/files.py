import os
import uuid

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)

def save_upload_bytes(upload_dir: str, filename: str, content: bytes) -> str:
    ensure_dir(upload_dir)
    ext = os.path.splitext(filename)[-1].lower()
    token = str(uuid.uuid4())
    safe_name = f"{token}{ext if ext else ''}"
    out_path = os.path.join(upload_dir, safe_name)
    with open(out_path, "wb") as f:
        f.write(content)
    return out_path

def public_url_for_file(base_url: str, static_mount: str, abs_or_rel_path: str) -> str:
    # We mount upload_dir at /static
    # Return a browser-loadable URL
    rel = abs_or_rel_path.replace("\\", "/")
    if rel.startswith(static_mount):
        return f"{base_url}{rel}"
    return f"{base_url}{static_mount}/{rel.split('/')[-1]}"