import os
import tarfile
import requests
import json
from tqdm import tqdm

def download_file(doi, filepath):
    url = 'https://api.datacite.org/dois/'+doi+'/media'
    r = requests.get(url).json()
    found = False
    for item in r['data']:
        url = item['attributes']['url']
        # Extract the last part of the URL after the last slash
        url_filename = url.split('/')[-1]
        if url_filename == filename:
            netcdf_url = url
            found = True
            break  # Exit the loop since we've found the matching filename
    if not found:
        print("Error: No matching filename found in the data.")
        return None
    r = requests.get(netcdf_url,stream=True)
    
    os.makedirs(data_path_root, exist_ok=True)
    full_path = os.path.join(data_path_root, filename)
    
    #Download file with progress bar
    if r.status_code == 403:
        print("File Unavailable")
        return None
    if 'content-length' not in r.headers:
        print("Did not get file")
        return None
    else:
        with open(full_path, 'wb') as f:
            total_length = int(r.headers.get('content-length'))
            for chunk in tqdm(r.iter_content(chunk_size=1024), total=total_length/1024, unit="KB"):
                if chunk:
                    f.write(chunk)
        return full_path
    
def download_and_extract(doi, filename, data_path_root, final_path):
    full_path = download_file(doi, filename, data_path_root)
    filepath_final = final_path
    if full_path:
        os.makedirs(final_path, exist_ok = True)
        if filename.endswith('.tar.gz'):
            try:
                with tarfile.open(full_path, 'r:gz') as tar:
                    tar.extractall(path=final_path)
                print(f"Extraction complete: {filename}. Contents are in {final_path}")
            except Exception as e:
                print(f"Error extracting the file: {e}")
        else:
            ### Insert code here to support other file formats
            filepath_final = full_path
    return filepath_final
