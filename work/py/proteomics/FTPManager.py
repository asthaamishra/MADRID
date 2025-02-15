"""
This file is responsible downloading data found at FTP links

TODO: Find a way to mark a file as "downloaded"
    - Keep a list of file names in a ".completd" hidden folder?
"""
from FileInformation import FileInformation
from FileInformation import clear_print
from ftplib import FTP
import multiprocessing
from multiprocessing.sharedctypes import Synchronized
import numpy as np
import time
from urllib.parse import urlparse


def ftp_client(host: str, max_attempts: int = 3) -> FTP:
    """
    This class is responsible for creating a "client" connection
    """
    connection_successful: bool = False
    attempt_num: int = 1
    # Attempt to connect, throw error if unable to do so
    while not connection_successful and attempt_num <= max_attempts:
        try:
            ftp_client: FTP = FTP(host, user="anonymous", passwd="guest")
            connection_successful = True
        except ConnectionResetError:
            
            # Make sure this print statement is on a new line on the first error
            if attempt_num == 1:
                print("")
            
            # Line clean: https://stackoverflow.com/a/5419488/13885200
            clear_print(f"Attempt {attempt_num} of {max_attempts} failed to connect")
            attempt_num += 1
            time.sleep(5)
    if not connection_successful:
        print("")
        raise ConnectionResetError("Could not connect to FTP server")
    
    return ftp_client


class Reader:
    def __init__(self, root_link: str, file_extensions: list[str]):
        """
        This class is responsible for reading data about root FTP links
        """
        self._root_link: str = root_link
        self._extensions: list[str] = file_extensions
        self._files: list[str] = []
        self._file_sizes: list[int] = []
        
        self._get_info()
        
    def _get_info(self):
        """
        This function is responsible for getting all files under the root_link
        """
        scheme: str = urlparse(self._root_link).scheme
        host = urlparse(self._root_link).hostname
        folder = urlparse(self._root_link).path
        
        client: FTP = ftp_client(host=host)
        
        for file_path in client.nlst(folder):
            if file_path.endswith(tuple(self._extensions)):
                download_url: str = f"{scheme}://{host}{file_path}"
                self._files.append(download_url)
                self._file_sizes.append(client.size(file_path))

    @property
    def files(self):
        return self._files
    
    @property
    def file_sizes(self):
        return self._file_sizes


class Download:
    def __init__(
        self,
        file_information: list[FileInformation],
        core_count: int = 1,
    ):
        """
        This function is responsible for downloading items from the FTP urls passed into it
        """
        self._file_information: list[FileInformation] = file_information
        self._core_count: int = min(core_count, 2)  # Limit to 2 downloads at a time
        self._download_counter: Synchronized = multiprocessing.Value("i", 1)

        # Find files to download
        self.download_data_wrapper()

    def download_data_wrapper(self):
        """
        This function is responsible for using multiprocessing to download as many files at once in parallel
        """
        # Calculate the number of cores ot use
        # We are going to use half the number of cores, OR the number of files to download, whichever is less
        # Otherwise if user specified the number of cores, use that
        print("Starting file download")

        # Split list into chunks to process separately
        file_chunks: list[FileInformation] = np.array_split(self._file_information, self._core_count)

        # Create a list of jobs
        jobs: list[multiprocessing.Process] = []
        
        for i, information in enumerate(file_chunks):

            # Append a job to the list
            job = multiprocessing.Process(
                target=self.download_data,
                args=(information,),
            )
            jobs.append(job)

        [job.start() for job in jobs]       # Start jobs
        [job.join() for job in jobs]        # Wait for jobs to finish
        [job.terminate() for job in jobs]   # Terminate jobs

    def download_data(self, file_information: list[FileInformation]):

        # Start processing the files
        for i, information in enumerate(file_information):
        
            # Define FTP-related variables
            parsed_url = urlparse(information.download_url)
            host = parsed_url.hostname
            path = parsed_url.path

            # Convert file_size from byte to MB
            size_mb: int = round(information.file_size / (1024**2))

            # Connect to the host, login, and download the file
            client: FTP = ftp_client(host=host)
            # client: FTP = FTP(host=host, user="anonymous", passwd="guest")

            # Get the lock and print file info
            self._download_counter.acquire()
            print(f"Started download {self._download_counter.value:02d} / {len(self._file_information):02d} ({size_mb} MB) - {information.raw_file_name}")  # fmt: skip
            self._download_counter.value += 1
            self._download_counter.release()

            # Download the file
            information.raw_file_path.parent.mkdir(parents=True, exist_ok=True)
            client.retrbinary(f"RETR {path}", open(information.raw_file_path, "wb").write)
            client.quit()


if __name__ == "__main__":
    print("Use the proteomics_preprocess.py file")
