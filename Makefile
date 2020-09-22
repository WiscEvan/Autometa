hello:
	@echo "Please inspect Makefile for list of commands"

clean:
	rm -rf htmlcov

cmd = --durations=0 --cov=autometa --emoji --cov-report html

test-wip: tests/test_data.json
	python -m pytest -m "wip" ${cmd}

test: tests/test_data.json
	python -m pytest ${cmd}

.Phony:
	test test-wip hello clean develop
