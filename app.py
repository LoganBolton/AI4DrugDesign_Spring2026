import gradio as gr
import os
from openai import OpenAI
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Initialize OpenAI client
client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))

def chat(message, history):
    """
    Chat function that sends messages to OpenAI API

    Args:
        message: The user's current message
        history: List of previous messages in the conversation

    Returns:
        The assistant's response
    """
    # Convert Gradio history format to OpenAI messages format
    messages = []
    for user_msg, assistant_msg in history:
        messages.append({"role": "user", "content": user_msg})
        messages.append({"role": "assistant", "content": assistant_msg})

    # Add the current message
    messages.append({"role": "user", "content": message})

    # Call OpenAI API
    response = client.chat.completions.create(
        model="gpt-5-nano-2025-08-07",
        messages=messages
    )

    return response.choices[0].message.content

# Create Gradio ChatInterface
demo = gr.ChatInterface(
    fn=chat,
    title="AI4DrugDesign Chat",
    description="Chat with GPT-5-nano-2025-08-07",
    examples=[
        "What is drug design?",
        "Explain molecular docking",
        "What are the challenges in AI-driven drug discovery?"
    ]
)

demo.launch()

