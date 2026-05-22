from bs4 import BeautifulSoup
import re

with open("sgm.htm", "r", encoding="utf-8") as f:
    soup = BeautifulSoup(f, "html.parser")

rows = soup.find_all("tr")
results = []

for row in rows:
    cells = row.find_all("td")
    if not cells:
        continue

    first_cell = cells[0].get_text(strip=True)
    match = re.match(r"\[(\d+)\]", first_cell)
    if not match:
        continue

    sg_number = int(match.group(1))
    names = [cell.get_text(strip=True) for cell in cells[1:] if cell.get_text(strip=True) and cell.get_text(strip=True) != ' ']

    for i, name in enumerate(names):
        results.append(str(sg_number + i))
        results.append(name)

# Write output with trailing newline
with open("spacegroup_names.txt", "w", encoding="utf-8") as out:
    out.write("\n".join(results) + "\n")

print(f"Extracted {len(results)//2} space group names to spacegroup_names.txt")

