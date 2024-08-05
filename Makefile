.PHONY: default
default:
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@echo "  protolint-fix  Run protolint with --fix option"
	@echo "  protolint-check  Run protolint to check the proto files"


.PHONY: protolint-fix
protolint-fix:
	protolint lint -fix .

.PHONY: protolint-check
protolint-check:
	protolint lint .
