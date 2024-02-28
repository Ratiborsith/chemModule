import requests
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
def extract_chemical_formulas(html_content):
    soup = BeautifulSoup(html_content, 'html.parser')
    chemical_formulas = []
    # Find all li elements with itemprop="match" attribute
    li_elements = soup.find_all('li', itemprop='match')
    for li in li_elements:
        # Find span element with itemprop="domain" and check if it contains "Chemical compound"
        domain_span = li.find('span', itemprop='domain')
        if domain_span and domain_span.text.strip() == "Chemical compound":
            # If it's a chemical compound, extract the name
            name_span = li.find('span', itemprop='name')
            if name_span:
                chemical_formulas.append(name_span.text.strip())
    return chemical_formulas

# Main function
def main():
    url = "https://patents.google.com/patent/TWI694074B"
    webpage_content = fetch_webpage(url)
    if webpage_content:
        chemical_formulas = extract_chemical_formulas(webpage_content)
        for formula in chemical_formulas:
            print("Chemical formula:", formula)

if __name__ == "__main__":
    main()
