#!/usr/bin/env python
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import sys


def read_urls(path):
    return Path(path).read_text().splitlines()


def main():
    urls = read_urls(sys.argv[1])

    download_dir = Path("downloads")
    download_dir.mkdir(exist_ok=True, parents=True)

    already_downloaded = []
    url2download = {}
    for url in urls:
        fname = url.rpartition("/")[-1]
        download_path = download_dir / fname
        if download_path.exists():
            already_downloaded.append(url)
        else:
            url2download[url] = download_path

    print(f"Already downloaded: \n{'\n'.join(already_downloaded)}\n")
    if any(url2download):
        print(f"Downloading: \n {'\n'.join(url2download.keys())}\n")
    else:
        print("Nothing to download")

    for url in url2download:
        print(f"Requesting {url}")
        req = Request(url)
        try:
            response = urlopen(req)
        except HTTPError as e:
            print("The server couldn't fulfill the request.")
            print("Error code: ", e.code)
        except URLError as e:
            print("We failed to reach a server.")
            print("Reason: ", e.reason)
        else:
            download_path = url2download[url]
            print(f"Saving {url} to {download_path}")
            with open(download_path, "wb") as f:
                f.write(response.read())


if __name__ == "__main__":
    main()
