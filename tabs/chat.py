import gradio as gr
import os
from openai import OpenAI

client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))


def _chat(message, history):
    messages = []
    for user_msg, assistant_msg in history:
        messages.append({"role": "user", "content": user_msg})
        messages.append({"role": "assistant", "content": assistant_msg})

    messages.append({"role": "user", "content": message})

    response = client.chat.completions.create(
        model="gpt-5-nano-2025-08-07",
        messages=messages,
    )

    return response.choices[0].message.content


def create_tab():
    with gr.Tab("Chat"):
        gr.ChatInterface(
            fn=_chat,
            title="AI4DrugDesign Chat",
            description="Chat with GPT-5-nano-2025-08-07",
            examples=[
                "What is drug design?",
                "Explain molecular docking",
                "What are the challenges in AI-driven drug discovery?",
            ],
        )
