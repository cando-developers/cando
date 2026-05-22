import os
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
        soup = BeautifulSoup(f, "html.parser")

    triplets = []
    for pre in soup.find_all("pre"):
        lines = pre.get_text().splitlines()
        for line in lines:
            line = line.strip()
            if line.count(",") == 2 and all(c in "xyzXYZ+-*/., 0123456789" for c in line.replace(",", "")):
                triplets.append(line)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(triplets))

    print(f"Extracted {len(triplets)} triplets for SG {sgnum}")
