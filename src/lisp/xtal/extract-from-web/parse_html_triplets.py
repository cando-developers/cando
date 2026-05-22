import os
import re
from bs4 import BeautifulSoup

input_dir = "symmetry_triplets"
output_dir = "parsed_triplets"
os.makedirs(output_dir, exist_ok=True)

for filename in os.listdir(input_dir):
    if not filename.endswith(".txt"):
        continue

    sgnum = filename.replace("sg", "").replace(".txt", "")
    input_path = os.path.join(input_dir, filename)
    output_path = os.path.join(output_dir, f"sg{sgnum}_triplets.txt")

    with open(input_path, "r", encoding="utf-8") as f:
        html = f.read()

    soup = BeautifulSoup(html, "html.parser")
    triplets = []

    for tag in soup.find_all("big"):
        match = re.search(r'\d+\s+(.+)', tag.text.strip())
        if match:
            triplet = match.group(1)
            if triplet.count(',') == 2:
                triplets.append(triplet)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(triplets))

    print(f"Extracted {len(triplets)} triplets for SG {sgnum}")
