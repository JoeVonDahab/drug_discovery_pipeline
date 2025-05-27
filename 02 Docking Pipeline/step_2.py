import json
from huggingface_hub import hf_hub_download
# pip install accelerate transformers
from transformers import AutoTokenizer, AutoModelForCausalLM
# pip install transformers
from transformers import pipeline
# pip install accelerate transformers
from transformers import AutoTokenizer, AutoModelForCausalLM
# pip install transformers
from transformers import pipeline

# Load prompt template for tasks from TDC
tdc_prompts_filepath = hf_hub_download(
    repo_id="google/txgemma-27b-chat",
    filename="tdc_prompts.json",
)
with open(tdc_prompts_filepath, "r") as f:
    tdc_prompts_json = json.load(f)

# Set example TDC task and input
task_name = "BBB_Martins"
input_type = "{Drug SMILES}"
drug_smiles = "CN1C(=O)CN=C(C2=CCCCC2)c2cc(Cl)ccc21"

# Construct prompt using template and input drug SMILES string
TDC_PROMPT = tdc_prompts_json[task_name].replace(input_type, drug_smiles)
print(TDC_PROMPT)



# Load model directly from Hugging Face Hub
tokenizer = AutoTokenizer.from_pretrained("google/txgemma-27b-chat")
model = AutoModelForCausalLM.from_pretrained(
    "google/txgemma-27b-chat",
    device_map="auto",
)

# Formatted TDC prompt (see "Formatting prompts for therapeutic tasks" section above)
prompt = TDC_PROMPT

# Prepare tokenized inputs
input_ids = tokenizer(prompt, return_tensors="pt").to("cuda")

# Generate response
outputs = model.generate(**input_ids, max_new_tokens=8)
print(tokenizer.decode(outputs[0], skip_special_tokens=True))



# Instantiate a text generation pipeline using the model
pipe = pipeline(
    "text-generation",
    model="google/txgemma-27b-chat",
    device="cuda",
)

# Formatted TDC prompt (see "Formatting prompts for therapeutic tasks" section above)
prompt = TDC_PROMPT

# Generate response
outputs = pipe(prompt, max_new_tokens=8)
response = outputs[0]["generated_text"]
print(response)



# Load model directly from Hugging Face Hub
tokenizer = AutoTokenizer.from_pretrained("google/txgemma-27b-chat")
model = AutoModelForCausalLM.from_pretrained(
    "google/txgemma-27b-chat",
    device_map="auto",
)

# Formatted TDC prompt (see "Formatting prompts for therapeutic tasks" section above)
prompt = TDC_PROMPT

# Format prompt in the conversational format
messages = [
    { "role": "user", "content": prompt}
]

# Apply the tokenizer's built-in chat template
chat_prompt = tokenizer.apply_chat_template(messages, tokenize=False, add_generation_prompt=True)

inputs = tokenizer.encode(chat_prompt, add_special_tokens=False, return_tensors="pt")
outputs = model.generate(input_ids=inputs.to("cuda"), max_new_tokens=8)
response = tokenizer.decode(outputs[0, len(inputs[0]):], skip_special_tokens=True)
print(response)

messages.extend([
    { "role": "assistant", "content": response },
    { "role": "user", "content": "Explain your reasoning based on the molecule structure." },
])

chat_prompt = tokenizer.apply_chat_template(messages, tokenize=False, add_generation_prompt=True)
inputs = tokenizer.encode(chat_prompt, add_special_tokens=False, return_tensors="pt")
outputs = model.generate(input_ids=inputs.to("cuda"), max_new_tokens=512)
response = tokenizer.decode(outputs[0, len(inputs[0]):], skip_special_tokens=True)
print(response)



# Instantiate a text generation pipeline using the model
pipe = pipeline(
    "text-generation",
    model="google/txgemma-27b-chat",
    device="cuda",
)

# Formatted TDC prompt (see "Formatting prompts for therapeutic tasks" section above)
prompt = TDC_PROMPT

# Format prompt in the conversational format for initial turn
messages = [
    { "role": "user", "content": prompt}
]

# Generate response for initial turn
outputs = pipe(messages, max_new_tokens=8)
print(outputs[0]["generated_text"][-1]["content"].strip())

# Append user prompt for an additional turn
messages = outputs[0]["generated_text"]
messages.append(
    { "role": "user", "content": "Explain your reasoning based on the molecule structure." }
)

# Generate response for additional turn
outputs = pipe(messages, max_new_tokens=512)
print(outputs[0]["generated_text"][-1]["content"].strip())
