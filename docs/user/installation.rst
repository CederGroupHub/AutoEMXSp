=========================
Installation
=========================

On a normal desktop computer, installation typically takes at most a few
minutes.

You can install the package using pip:

.. code-block:: bash

   pip install autoemx

To use AutoEMX through its GUI (acquisition, quantification, analysis and single-spectrum
fits), create the double-click app once:

.. code-block:: bash

   python -m autoemx.gui --create-launcher

See :ref:`gui_getting_started_tutorial`.

Optional particle segmentation models (see :ref:`particle_size_statistics`), which
segment touching particles separately:

.. code-block:: bash

   pip install "autoemx[segmentation]"                          # SAM2
   python -m autoemx.utils.install_segmentation micro-sam      # micro-sam

pip then warns that micro-sam requires PyQt6, segment-anything-py and trackastra: ignore
it, they are left out on purpose.

With an editable installation (see below), install SAM2 from the repository folder
instead: ``pip install -e ".[segmentation]"``.


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
