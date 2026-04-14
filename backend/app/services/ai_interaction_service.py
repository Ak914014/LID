import requests
import json

OLLAMA_URL = "http://localhost:11434/api/generate"
MODEL = "mistral"


def call_ollama(prompt: str):
    try:
        res = requests.post(
            OLLAMA_URL,
            json={
                "model": "mistral",
                "prompt": prompt,
                "stream": False
            },
            timeout=20
        )
        return res.json().get("response", "")
    except Exception:
        return None


def rank_and_filter_interactions(interactions):
    """
    interactions: List[dict]
    Expected keys:
    - type
    - distance
    - residue
    - ligand_atom
    """

    # ---- Rule-based pre-filter (fast + safe) ----
    allowed = {"hbond", "hydrophobic", "pi_stack", "salt_bridge"}

    filtered = [
        i for i in interactions
        if i.get("type") in allowed and i.get("distance", 10) < 4.5
    ]

    if len(filtered) <= 8:
        return filtered

    # ---- AI Ranking ----
    prompt = f"""
    You are an expert in protein-ligand binding.

    Given these interactions:
    {json.dumps(filtered, indent=2)}

    Select the TOP 8 most important interactions.

    Prefer:
    - hydrogen bonds
    - pi stacking
    - salt bridges
    - shorter distances

    Return ONLY JSON list.
    """

    response = call_ollama(prompt)

    try:
        ranked = json.loads(response)
        return ranked[:8]
    except Exception:
        return filtered[:8]