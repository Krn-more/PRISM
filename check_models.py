import os
import google.generativeai as genai

key_path = os.path.join(os.path.dirname(__file__), "..", "api_key.txt")
if os.path.exists(key_path):
    with open(key_path, 'r') as f:
        api_key = f.read().strip()
else:
    api_key = os.environ.get("GEMINI_API_KEY", "")

if not api_key:
    print("API key not found")
else:
    genai.configure(api_key=api_key)
    try:
        for m in genai.list_models():
            if 'generateContent' in m.supported_generation_methods:
                print(m.name)
    except Exception as e:
        print(f"Error listing models: {e}")
