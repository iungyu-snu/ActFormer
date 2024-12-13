import requests
from typing import Dict, List
import json
import sys
import xml.etree.ElementTree as ET
import traceback

class SabioRKAnalyzer:
    def __init__(self):
        self.base_url = "http://sabiork.h-its.org/sabioRestWebServices"
        self.uniprot_url = "https://rest.uniprot.org/uniprotkb"
        
    def get_kinetics_by_ec(self, ec_number: str) -> List[Dict]:
        """Query SABIO-RK for kinetic parameters matching an EC number"""
        # First, search for kinetic law IDs
        search_endpoint = f"{self.base_url}/searchKineticLaws/sbml"
        
        # Parameters for the search query
        search_params = {
            "q": f"ECNumber:{ec_number}",
            "output": "xml"
        }
        
        print(f"\nDebug: Searching kinetic laws")
        print(f"Debug: Search endpoint: {search_endpoint}")
        print(f"Debug: Search parameters: {search_params}")
        
        search_response = requests.get(search_endpoint, params=search_params)
        print(f"Debug: Search response status: {search_response.status_code}")
        
        if search_response.status_code != 200:
            raise Exception(f"Search failed: {search_response.status_code}")
            
        # Extract kinetic law IDs from search results
        try:
            root = ET.fromstring(search_response.text)
            kinlaw_ids = []
            
            # Find all kineticLawID elements
            for reaction in root.findall('.//*'):
                if reaction.tag.endswith('kineticLawID'):
                    kid = reaction.text
                    print(f"Debug: Found kinetic law ID: {kid}")
                    kinlaw_ids.append(kid)
            
            print(f"\nDebug: Found {len(kinlaw_ids)} kinetic law IDs")
            
            if not kinlaw_ids:
                print(f"No kinetic laws found for EC {ec_number}")
                return []
                
            # Now get detailed data for each kinetic law ID
            results = []
            for kid in kinlaw_ids:
                print(f"\nDebug: Fetching details for kinetic law ID: {kid}")
                detail_endpoint = f"{self.base_url}/kineticLaws"
                detail_params = {
                    "kinlawids": kid,
                    "output": "xml"
                }
                detail_response = requests.get(detail_endpoint, params=detail_params)
                print(f"Debug: Detail response status: {detail_response.status_code}")
                print(f"Debug: Detail response content preview: {detail_response.text[:500]}...")
                
                if detail_response.status_code == 200:
                    try:
                        detail_root = ET.fromstring(detail_response.text)
                        # Print out the XML structure
                        print("Debug: XML structure:")
                        self._print_xml_structure(detail_root)
                        
                        entry_data = self._parse_kinetic_entry(detail_root)
                        if entry_data:
                            print(f"Debug: Successfully parsed entry data")
                            results.append(entry_data)
                        else:
                            print(f"Debug: No valid data found in entry")
                    except Exception as e:
                        print(f"Debug: Error parsing entry {kid}: {str(e)}")
                else:
                    print(f"Debug: Failed to fetch details for ID {kid}: {detail_response.status_code}")
                    
            return results
            
        except Exception as e:
            print(f"Debug: Error processing search results: {str(e)}")
            raise
            
    def _print_xml_structure(self, element, level=0):
        """Helper method to print XML structure"""
        print("  " * level + f"Tag: {element.tag}")
        if element.attrib:
            print("  " * level + f"Attributes: {element.attrib}")
        if element.text and element.text.strip():
            print("  " * level + f"Text: {element.text.strip()[:100]}")
        for child in element:
            self._print_xml_structure(child, level + 1)
            
    def get_protein_sequence(self, uniprot_id: str) -> str:
        """Fetch protein sequence from UniProt"""
        try:
            response = requests.get(f"{self.uniprot_url}/{uniprot_id}.fasta")
            if response.status_code == 200:
                # Skip the header line and join the sequence lines
                sequence = ''.join(response.text.split('\n')[1:])
                print(f"Debug: Found sequence of length {len(sequence)} for {uniprot_id}")
                return sequence
            else:
                print(f"Debug: Failed to fetch sequence for {uniprot_id}")
                return None
        except Exception as e:
            print(f"Debug: Error fetching sequence: {str(e)}")
            return None

    def parse_mutation(self, protein_name: str) -> List[tuple]:
        """Parse mutation information from protein name"""
        mutations = []
        if 'mutant' in protein_name.lower():
            # Extract mutation codes (e.g., Q240E)
            parts = protein_name.split()
            for part in parts:
                # Match pattern like Q240E
                if len(part) >= 3 and part[0].isalpha() and part[-1].isalpha() and part[1:-1].isdigit():
                    wt = part[0]  # Wild type amino acid
                    pos = int(part[1:-1])  # Position
                    mut = part[-1]  # Mutant amino acid
                    mutations.append((wt, pos, mut))
        return mutations

    def verify_and_modify_sequence(self, sequence: str, mutations: List[tuple]) -> str:
        """Verify mutations against sequence and modify if valid"""
        if not mutations:
            return sequence
        
        modified_seq = list(sequence)
        errors = []
        
        for wt, pos, mut in mutations:
            # Convert position to 0-based index
            idx = pos - 1
            
            # Verify position is within sequence
            if idx >= len(sequence):
                errors.append(f"Error: Mutation position {pos} is outside sequence length {len(sequence)}")
                continue
            
            # Verify wild type amino acid
            if sequence[idx] != wt:
                errors.append(f"Error: Expected {wt} at position {pos} but found {sequence[idx]}")
                continue
            
            # Apply mutation
            modified_seq[idx] = mut
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return "".join(modified_seq)

    def _parse_kinetic_entry(self, root: ET.Element) -> Dict:
        """Parse a single kinetic law entry"""
        try:
            # Define namespaces
            ns = {
                'sbml': 'http://www.sbml.org/sbml/level3/version1/core',
                'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
                'sabiork': 'http://sabiork.h-its.org',
                'bqbiol': 'http://biomodels.net/biology-qualifiers/'
            }
            
            entry = {}
            
            # Get enzyme information
            for species in root.findall('.//sbml:species', ns):
                name = species.get('name', '')
                if 'enzyme' in name.lower() or 'catalyst' in name.lower():
                    entry['protein_name'] = name
                    # Look for UniProt ID
                    for li in species.findall('.//rdf:li', ns):
                        resource = li.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource', '')
                        if 'uniprot' in resource.lower():
                            entry['uniprot_id'] = resource.split('/')[-1]
                            print(f"Debug: Found UniProt ID: {entry['uniprot_id']}")
                    break
            
            # Get substrate information
            substrates = []
            for species in root.findall('.//sbml:species', ns):
                name = species.get('name', '')
                if name and not ('enzyme' in name.lower() or 'catalyst' in name.lower()):
                    substrates.append(name)
                    print(f"Debug: Found substrate: {name}")
            entry['substrate'] = '; '.join(filter(None, substrates))
            print(f"Debug: Combined substrates: {entry['substrate']}")
            
            # Get kinetic parameters from kineticLaw
            for kinetic_law in root.findall('.//sbml:kineticLaw', ns):
                # Look for Vmax and Km parameters
                for param in kinetic_law.findall('.//sbml:localParameter', ns):
                    param_id = param.get('id', '').lower()
                    param_name = param.get('name', '').lower()
                    param_value = param.get('value')
                    param_units = param.get('units', 'N/A')
                    
                    print(f"Debug: Found parameter: {param_id} ({param_name}) = {param_value} {param_units}")
                    
                    if 'vmax' in param_id or 'vmax' in param_name:
                        entry['vmax'] = param_value
                        entry['vmax_unit'] = param_units
                    elif 'km' in param_id or 'km' in param_name:
                        entry['km'] = param_value
                        entry['km_unit'] = param_units
                
                # Get experimental conditions
                conditions = kinetic_law.find('.//sabiork:experimentalConditions', ns)
                if conditions is not None:
                    temp = conditions.find('.//sabiork:startValueTemperature', ns)
                    ph = conditions.find('.//sabiork:startValuepH', ns)
                    buffer = conditions.find('.//sabiork:buffer', ns)
                    
                    if temp is not None:
                        temp_unit = conditions.find('.//sabiork:temperatureUnit', ns)
                        entry['temperature'] = f"{temp.text} {temp_unit.text if temp_unit is not None else '°C'}"
                    if ph is not None:
                        entry['ph'] = ph.text
                    if buffer is not None:
                        entry['buffer'] = buffer.text
                    
                # Calculate kcat/KM if we have both values
                if 'vmax' in entry and 'km' in entry:
                    try:
                        vmax = float(entry['vmax'])
                        km = float(entry['km'])
                        entry['kcat_km'] = str(vmax / km)
                        entry['unit'] = f"{entry['vmax_unit']}/{entry['km_unit']}"
                        print(f"Debug: Calculated kcat/KM: {entry['kcat_km']} {entry['unit']}")
                    except (ValueError, ZeroDivisionError) as e:
                        print(f"Debug: Error calculating kcat/KM: {str(e)}")
            
            if 'uniprot_id' in entry:
                # Fetch protein sequence
                sequence = self.get_protein_sequence(entry['uniprot_id'])
                if sequence:
                    # Parse mutations from protein name
                    mutations = self.parse_mutation(entry['protein_name'])
                    try:
                        # Verify and modify sequence if there are mutations
                        modified_sequence = self.verify_and_modify_sequence(sequence, mutations)
                        entry['sequence'] = modified_sequence
                        if mutations:
                            entry['mutations'] = [f"{wt}{pos}{mut}" for wt, pos, mut in mutations]
                    except ValueError as e:
                        entry['sequence'] = sequence
                        entry['sequence_error'] = str(e)
                else:
                    print(f"Debug: Could not fetch sequence for {entry['uniprot_id']}")
                    return None
            
            # Only return entry if it has required fields
            required_fields = ['sequence', 'substrate', 'kcat_km', 'unit']
            missing_fields = [f for f in required_fields if f not in entry]
            if missing_fields:
                print(f"Debug: Missing required fields: {missing_fields}")
                return None
            return entry
            
        except Exception as e:
            print(f"Debug: Error parsing kinetic entry: {str(e)}")
            traceback.print_exc()
            return None

    def save_results(self, ec_number: str, entries: List[Dict], filename: str = "result.txt"):
        """Save kinetic parameters to a file"""
        with open(filename, 'w') as f:
            f.write(f"Kinetic Parameters for EC {ec_number}:\n")
            f.write("--------------------------------------------------\n\n")
            
            if not entries:
                f.write("No kinetic data found\n")
                return
            
            for i, entry in enumerate(entries, 1):
                f.write(f"Entry {i}:\n")
                f.write(f"Protein Name: {entry.get('protein_name', 'Unknown')}\n")
                if 'uniprot_id' in entry:
                    f.write(f"UniProt ID: {entry['uniprot_id']}\n")
                if 'mutations' in entry:
                    f.write(f"Mutations: {', '.join(entry['mutations'])}\n")
                if 'sequence_error' in entry:
                    f.write(f"Sequence Error: {entry['sequence_error']}\n")
                if 'sequence' in entry:
                    f.write(f"Sequence: {entry['sequence']}\n")
                f.write(f"Substrate(s): {entry['substrate']}\n")
                f.write(f"kcat/KM: {entry['kcat_km']} {entry['unit']}\n")
                if 'temperature' in entry:
                    f.write(f"Temperature: {entry['temperature']}\n")
                if 'ph' in entry:
                    f.write(f"pH: {entry['ph']}\n")
                if 'buffer' in entry:
                    f.write(f"Buffer: {entry['buffer']}\n")
                f.write("\n")

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 data.py <EC_number>")
        sys.exit(1)
        
    analyzer = SabioRKAnalyzer()
    ec_number = sys.argv[1]
    entries = analyzer.get_kinetics_by_ec(ec_number)
    
    # Print to console
    print("\nKinetic Parameters for EC {}:".format(ec_number))
    print("--------------------------------------------------\n")
    if entries:
        print(f"Total entries shown: {len(entries)}")
        for entry in entries:
            print(json.dumps(entry, indent=2))
    else:
        print("No kinetic data found for EC {}".format(ec_number))
    
    # Save to file
    analyzer.save_results(ec_number, entries)
    print("\nResults have been saved to result.txt")

if __name__ == "__main__":
    main()
