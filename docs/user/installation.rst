=========================
Installation
=========================

On a normal desktop computer, installation typically takes at most a few
minutes.

You can install the package using pip:

.. code-block:: bash

   pip install autoemx


Alternatively, to ensure you are constantly up to date with the package during minor developments, you can clone the GitHub repository and use editable installation. This allows you to access the latest features or bug fixes that are not yet included in the latest release version.
To do this, follow these steps:

1. **Navigate to the Directory:**
   Change into the directory you want the package to be cloned:

.. code-block:: bash
   
    cd your/path/

2. **Clone the Repository:**
   Create a local copy of the repository at the current directory:
   
.. code-block:: bash

    git clone https://github.com/CederGroupHub/AutoEMX.git

3. **Install the Package in Editable Mode:**

.. code-block:: bash
   
    pip install -e .


To ensure that you always have the latest version of the package, complete with all new features and fixes, periodically run:

1. **Navigate to the Directory:**

.. code-block:: bash
   
    cd your/path/to/AutoEMX
    
    
2. **Pull Updates:**
    Fetch latest package updates from the main branch of the repository:
    
.. code-block:: bash
   
    git pull origin main
