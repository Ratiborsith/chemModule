import requests
import re
from bs4 import BeautifulSoup

# Function to fetch webpage content
def fetch_webpage(url):
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        print("Failed to fetch webpage")
        return None

# Function to extract chemical formulas from HTML content
def extract_chemical_info(html_content):
    soup = BeautifulSoup(html_content, 'html.parser')
    chemical_info = {}
    # Find all li elements with itemprop="match" attribute
    li_elements = soup.find_all('li', itemprop='match')
    for li in li_elements:
        # Find span elements with relevant itemprop values
        domain_span = li.find('span', itemprop='domain')
        name_span = li.find('span', itemprop='name')
        smiles_span = li.find('span', itemprop='smiles')
        # Check if all relevant information is present and domain is "Chemical compound"
        if domain_span and domain_span.text.strip() == "Chemical compound" and name_span and smiles_span:
            domain = domain_span.text.strip()
            name = name_span.text.strip()
            smiles = smiles_span.text.strip()
            # Use regular expressions to identify types of formulas
            if re.match(r"^[^=;]*$", name):  # UIPAC
                chemical_info.setdefault("UIPAC", []).append(name)
            elif re.match(r'^[A-Za-z0-9@+\-\[\]\(\)\\\/%=#]+$', name):  # SMILES
                chemical_info.setdefault("SMILES", []).append(smiles)
            else:  # compoundName
                chemical_info.setdefault("compoundName", []).append(name)
    return chemical_info

# Main function
def main():
    url = "https://patents.google.com/patent/US8980879B2"
    webpage_content = fetch_webpage(url)
    if webpage_content:
        chemical_info = extract_chemical_info(webpage_content)
        # Print by the specified order
        i = 0
        for type_name in ["SMILES", "UIPAC", "compoundName"]:

            if type_name == "UIPAC":
                print(type_name + ":")
                for formula in chemical_info[type_name]:
                    print(formula)
                    i = i + 1


        print(f"Количество найденных формул: {i}")

if __name__ == "__main__":
    main()
