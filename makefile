install:
	pip3 uninstall -y karstnet
	pip3 install .

pep8:
	autopep8 --in-place karstnet/*.py
	autopep8 --in-place tests/*.py

test:
	pytest tests/*.py
	pycodestyle karstnet/*.py
	pycodestyle tests/*.py
