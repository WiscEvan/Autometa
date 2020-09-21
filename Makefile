
hello:
	echo "Please inspect Makefile for list of commands"


test-wip: tests/test_data.json
	python -m pytest -m "wip" --durations=0 --cov=autometa --emoji --cov-report html

test: tests/test_data.json
	python -m pytest --durations=0 --cov=autometa --emoji --cov-report html

.Phony:
	test test-wip hello
