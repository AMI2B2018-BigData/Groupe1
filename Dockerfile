# Use an official Python runtime as a parent image
FROM python:3-slim

# Copy the current directory contents into the container at /Projet_python_big_data
COPY . /Projet_python_big_data

# Set the working directory to /Projet_python_big_data
WORKDIR /Projet_python_big_data

# Run to build the image
RUN pip install --trusted-host pypi.python.org -r requirements.txt
RUN mkdir results

# Run Confo_global.py when the container launches
CMD ["python3", "Confo_global.py", "-ref", "pab21_structure_de_ref.pdb", "-conf", "pab21_prod_solute_500frames.pdb", "-domaines", "A1,A2,A3,A4,B"]
