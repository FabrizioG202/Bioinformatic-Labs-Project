import aiohttp
import asyncio
from Bio import SeqIO
from io import StringIO
import gzip
import shutil
import os
import urllib
 
async def get_pdb_resolution(pdb_id: str) -> float:
    """Get the resolution of a PDB structure"""
    # Set the URL for the Data API
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    async with aiohttp.ClientSession() as session:
        async with session.get(url) as response:
            # Check if the request was successful
            if response.status != 200:
                return None
            
            # Get the JSON data from the response
            data = await response.json()
            
            # Get the resolution of the structure
            return data.get('rcsb_entry_info', {}).get('resolution_combined', [None])[0]
    
    
# calls the get_pdb_resolution function ona list of PDB IDs
# and returns a dictionary mapping PDB IDs to resolution
async def get_pdb_resolution_list(pdb_ids: list) -> dict:
    """Get the resolution of a list of PDB structures"""
    tasks = [ get_pdb_resolution(pdb_id) for pdb_id in pdb_ids ]
    resolutions = await asyncio.gather(*tasks)
    
    return dict(zip(pdb_ids, resolutions))

# gets the sequence of a protein given its PDB ID
async def get_sequence(pdb_id: str, chain : str | None = None) -> str:
    """Get the sequence of a protein given its PDB ID"""
    # Set the URL for the Data API
    url = None
    if chain is not None:
        url = f'https://www.rcsb.org/fasta/chain/{pdb_id}.{chain}/download'
    else:  
        url = f'https://www.rcsb.org/fasta/entry/{pdb_id}/download'

    async with aiohttp.ClientSession() as session:
        async with session.get(url) as response:
            # Check if the request was successful
            if response.status != 200:
                print(f"Failed to get sequence for {pdb_id}")
                return None
            
            # Get the sequence from the response
            data= await response.text()
            
            if chain is None:
                return data
            
            # parse the sequence and return the first one
            for record in SeqIO.parse(StringIO(data), "fasta"):
                return str(record.seq)
    
            return None
        
    
# calls the get_sequence function on a list of PDB IDs
# and returns a dictionary mapping PDB IDs to sequences
# takes in a list of string of the type pdb_id:chain, so they must  be formatted as such
async def get_sequence_list(pdb_and_chains: list[str]) -> dict:
    """Get the sequence of a list of proteins given their PDB IDs"""
    tasks = []
    for pdb_and_chain in pdb_and_chains:
        values = pdb_and_chain.split(':')
        if len(values) != 2:
            print(f"Invalid PDB ID and chain: {pdb_and_chain}")
            continue
        
        pdb_id, chain = values
        tasks.append(get_sequence(pdb_id, chain))
    sequences = await asyncio.gather(*tasks)
    
    return dict(zip(pdb_and_chains, sequences))


# Gets the sequence of a protein given its PDB ID and its chain id
# it opens the file and finds pdbid_chain in the file then returns the sequence (minus the header)
# be ware that the file is very large and this function is very slow 
# and that sequences can be more than one line long
# it takes in a list of pdb and chain ids like this: [pdbid:chain]
def get_sequences_from_file(pdbids : list[str], res_seq_path : str) -> str:
    """Get the sequence of a protein given its PDB ID"""
    
    clean_ids = []

    # convert anu id in the list to the format pdbid_chainid if it is not already in that format
    # the other format it might be in is pdbid:chainid
    for i in range(len(pdbids)):
        if len(pdbids[i].split(":")) == 2:
            clean_ids.append(pdbids[i].replace(":", "_"))
        else:
            clean_ids.append(pdbids[i])
    
    pdbids = dict()
    
    # open the file
    with open(res_seq_path, "r") as f:
        current_pdbid = None
        for line in (f):
            # check if the line is a header and if it (lowered) starts with pdbid_chainid
            if line.startswith(">"):
                for pdbid in clean_ids:
                    if line.lower().startswith(f">{pdbid.lower()}"):
                        current_pdbid = pdbid
                        break
                else:
                    current_pdbid = None                    
            else:
                # if the current pdbid is not none then we add the line to the current sequence
                if current_pdbid is not None:
                    if current_pdbid not in pdbids:
                        pdbids[current_pdbid] = ""
                    pdbids[current_pdbid] += line.strip()
    return pdbids


# download the pdb sequences from the pdb folder and saves it to a file called seq_res.fasta
# in the cache subdirectory (if it doesn't exist it will be created)
# the file is downloaded from https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
# and extracted to the cache subdirectory
# if name is provided it will be saved as name.fasta
# if it's not provided it will be saved as seq_res.fasta
# it returns the path to the file
# if overwrite is true it will overwrite the file if it already exists
# otherwise it will return the path to the file if it already exists
CACHE_FOLDER_NAME ="./cache"
def download_pdb_sequences(file_name : str = "seq_res", overwrite : bool = False) -> str:
    """Download the PDB sequences from the PDB folder and saves it to a file called seq_res.fasta"""
    

    # if the cache directory doesn't exist then create it
    if not os.path.exists(CACHE_FOLDER_NAME):
        os.mkdir(CACHE_FOLDER_NAME)
        
    final_file_name = f"{CACHE_FOLDER_NAME}/{file_name}.fasta"
    
    # if the file already exists and overwrite is false then return the path to the file
    if os.path.exists(final_file_name) and not overwrite:
        return final_file_name
    
    # download the file from the url
    print("Downloading PDB sequences...")
    url = "https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz"
    gzipped_file_name = f"{CACHE_FOLDER_NAME}/{file_name}.txt.gz"
    
    # download the file
    urllib.request.urlretrieve(url, gzipped_file_name)
    
    # extract the file
    with gzip.open(gzipped_file_name, 'rb') as f_in:
        with open(final_file_name, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    # delete the gzipped file
    os.remove(gzipped_file_name)
    
    return final_file_name
    