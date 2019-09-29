#################################################################################
# GLOBALS                                                                       #
#################################################################################
.DEFAULT_GOAL := check
PROJECT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
PROFILE = default
PROJECT_NAME = mlp
PYTHON_INTERPRETER = python3

################################################################################
# COMMANDS                                                                     #
# TO SHOW OUTPUT USE LOGGER=stdout                                             #
################################################################################


main:
	@echo "---> Starting Code to Train, Test and Validate"
	@$(PYTHON_INTERPRETER) src/api/main.py


setup: check_environment
	@echo "---> Running setup.."
	@conda env create -f environment.yml
	@cp -n .env.example .env
	@echo "---> To complete setup please run \n---> source activate cerberus"


install:
	@echo "---> Installing dependencies"
	@conda env update -f environment.yml


clean:
	find . -type f -name "*.py[co]" -delete
	find . -type d -name "__pycache__" -delete


lint:
	flake8 src


check_environment:
	@echo "---> Checking environment.."
	$(PYTHON_INTERPRETER) test_environment.py


autocorrect:
	@echo "---> Processing autocorrect"
	@autopep8 --in-place --aggressive --aggressive --global-config .flake8 $(shell find . -name '*.py')


console:
	@$(PYTHON_INTERPRETER)
