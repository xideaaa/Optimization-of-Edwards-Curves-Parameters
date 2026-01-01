CC       := clang
BASEC    := -O3 -Wall -Wextra -std=c11
CFLAGS   := $(BASEC)
LDFLAGS  :=

ifdef cost
CFLAGS   += -DFF_COUNT_OPS
endif

SRC_DIR     := src
EX_DIR      := examples
TST_DIR     := tests
BLD_DIR     := build
OBJ_DIR     := $(BLD_DIR)/obj
DATA_DIR    := data
REPORT_DIR  := reports
SCRIPTS_DIR := scripts

COMMON_SRC   := $(wildcard $(SRC_DIR)/*.c)
COMMON_SRC := $(filter-out \
    $(SRC_DIR)/curve_gen.c \
    $(SRC_DIR)/curve_select.c \
    $(SRC_DIR)/benchmark.c \
    $(SRC_DIR)/benchmark_iso.c, \
    $(COMMON_SRC))
COMMON_OBJ   := $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/src/%.o,$(COMMON_SRC))

EX_SRC       := $(wildcard $(EX_DIR)/*.c)
TST_SRC      := $(wildcard $(TST_DIR)/*.c)

EX_OBJ       := $(patsubst $(EX_DIR)/%.c,$(OBJ_DIR)/examples/%.o,$(EX_SRC))
TST_OBJ      := $(patsubst $(TST_DIR)/%.c,$(OBJ_DIR)/tests/%.o,$(TST_SRC))

EXE_EX       := $(patsubst $(EX_DIR)/%.c,$(BLD_DIR)/%,$(EX_SRC))
EXE_TST      := $(patsubst $(TST_DIR)/%.c,$(BLD_DIR)/%,$(TST_SRC))
ALL_BIN      := $(EXE_EX) $(EXE_TST)

GEN_BIN         := $(BLD_DIR)/curve_gen         # src/curve_gen.c
SELECT_BIN      := $(BLD_DIR)/curve_select      # src/curve_select.c
BENCH_BIN       := $(BLD_DIR)/ed25519_acc_bench # examples/ed25519_acc_bench.c
BENCH_CURVES_BIN:= $(BLD_DIR)/benchmark         # src/benchmark.c

CAND_CSV      := $(DATA_DIR)/candidates.csv
COUNTED_CSV   := $(DATA_DIR)/counted.csv
ACCEPTED_CSV  := $(DATA_DIR)/accepted.csv

SAGE          := sage
SAGE_SCRIPT   := $(SCRIPTS_DIR)/count_with_sage.sage

# Defaults for Sage
SMALL_BOUND ?= 100000
PRP ?= 1
SKIP_TWIST ?= 0
ALGO ?= sea

.PHONY: all
all: $(ALL_BIN) $(GEN_BIN) $(SELECT_BIN) $(BENCH_CURVES_BIN)

# examples/
$(EXE_EX): $(BLD_DIR)/%: $(OBJ_DIR)/examples/%.o $(COMMON_OBJ) | $(BLD_DIR)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

# tests/
$(EXE_TST): $(BLD_DIR)/%: $(OBJ_DIR)/tests/%.o $(COMMON_OBJ) | $(BLD_DIR)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

# main in src/ ---
$(GEN_BIN): $(OBJ_DIR)/src/curve_gen.o $(COMMON_OBJ) | $(BLD_DIR)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

$(SELECT_BIN): $(OBJ_DIR)/src/curve_select.o $(COMMON_OBJ) | $(BLD_DIR)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

$(BENCH_CURVES_BIN): $(OBJ_DIR)/src/benchmark.o $(COMMON_OBJ) | $(BLD_DIR)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

# compile
$(OBJ_DIR)/src/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)/src
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

$(OBJ_DIR)/examples/%.o: $(EX_DIR)/%.c | $(OBJ_DIR)/examples
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

$(OBJ_DIR)/tests/%.o: $(TST_DIR)/%.c | $(OBJ_DIR)/tests
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@


-include $(COMMON_OBJ:.o=.d) $(EX_OBJ:.o=.d) $(TST_OBJ:.o=.d) \
          $(OBJ_DIR)/src/curve_gen.d $(OBJ_DIR)/src/curve_select.d \
          $(OBJ_DIR)/src/benchmark.d


.PHONY: run-%
run-%: $(BLD_DIR)/%
	@echo "==> Running $<"
	@$<

.PHONY: clean
clean:
	rm -rf $(BLD_DIR)

$(BLD_DIR) \
$(OBJ_DIR) \
$(OBJ_DIR)/src \
$(OBJ_DIR)/examples \
$(OBJ_DIR)/tests \
$(DATA_DIR) \
$(REPORT_DIR):
	@mkdir -p $@


# STAGE 1 - candidate generator
# make gen N=1000 P=<hex> A=<hex> OUT=<cand.csv> (COMPLETE=<0/1>)
N ?= 100
P ?= 0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed   # Ed25519
A ?= -1
COMPLETE ?= 1
OUT ?= $(CAND_CSV)

.PHONY: gen
gen: $(DATA_DIR) $(GEN_BIN)
	@echo "==> Generating $(N) candidates into $(OUT)"
	@$(GEN_BIN) --p $(P) --a $(A) --count $(N) --complete $(COMPLETE) --out $(OUT)

# STAGE 2 (flexible paths)
# make count-file IN=<cand.csv> OUT=<counted.csv>
.PHONY: count-file
count-file:
ifneq ("$(wildcard $(SAGE_SCRIPT))","")
	@[ -n "$(IN)" ] && [ -n "$(OUT)" ] || (echo "Usage: make count-file IN=<cand.csv> OUT=<counted.csv>"; false)
	@mkdir -p $(dir $(OUT))
	@echo "==> Counting with Sage: IN=$(IN) -> OUT=$(OUT)"
	@SMALL_BOUND=$(SMALL_BOUND) PRP=$(PRP) SKIP_TWIST=$(SKIP_TWIST) ALGO=$(ALGO) \
	  $(SAGE) $(SAGE_SCRIPT) $(IN) $(OUT)
	@echo "==> Wrote $(OUT)"
else
	@echo "ERROR: Missing $(SAGE_SCRIPT). Please add it."
	@false
endif


# STAGE 3 - C selector
# make select-file IN=<counted.csv> OUT=<accepted.csv>
.PHONY: select-file
select-file: $(SELECT_BIN)
	@[ -n "$(IN)" ] && [ -n "$(OUT)" ] || (echo "Usage: make select-file IN=<counted.csv> OUT=<accepted.csv>"; false)
	@mkdir -p $(dir $(OUT))
	@echo "==> Selecting (C): IN=$(IN) -> OUT=$(OUT)"
	@$(SELECT_BIN) --in $(IN) --out $(OUT) --verbose 1



# STAGE 4 - benchmark (arthm) accepted.csv
# make bench-curves COST=<0/1> IN=<accepted.csv> OUT=<bench_results.csv>  ITERS=<no.>
ITERS ?= 20

.PHONY: bench-curves
bench-curves:
	@[ -n "$(IN)" ] && [ -n "$(OUT)" ] || (echo "Usage: make bench-curves IN=<accepted.csv> OUT=<bench.csv> [ITERS=N]"; false)
	@mkdir -p $(dir $(OUT))
	@echo "==> Benchmarking curves: IN=$(IN) -> OUT=$(OUT) ITERS=$(ITERS)"
	@$(MAKE) clean
	@$(MAKE) cost=1 $(BLD_DIR)/benchmark
	@$(BLD_DIR)/benchmark --in "$(IN)" --out "$(OUT)" --iters "$(ITERS)"


BENCH_ISO_CHAIN_BIN := $(BLD_DIR)/benchmark_iso

$(BENCH_ISO_CHAIN_BIN): $(OBJ_DIR)/src/benchmark_iso.o $(COMMON_OBJ) | $(BLD_DIR)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@


# STAGE 5 – benchmark (isogeny) accepted.csv
# make bench-iso COST=<0/1> IN=<accepted.csv> OUT=<bench_iso_results.csv>  ITERS=<no.>
.PHONY: bench-iso
bench-iso:
	@[ -n "$(IN)" ] && [ -n "$(OUT)" ] || \
	  (echo "Usage: make bench-iso IN=accepted.csv OUT=bench.csv ITERS=N"; false)
	@$(MAKE) clean
	@$(MAKE) cost=1 $(BENCH_ISO_CHAIN_BIN)
	@$(BENCH_ISO_CHAIN_BIN) --in "$(IN)" --out "$(OUT)" --iters "$(ITERS)"



.PHONY: help
help:
	@echo "Targets:"
	@echo "  make"
	@echo "  make gen N=1000 P=<hex> A=<hex> [COMPLETE=1] [OUT=path]"
	@echo "  make count-file IN=... OUT=... [SMALL_BOUND=2 PRP=1 SKIP_TWIST=0 ALGO=sea]"
	@echo "  make select-file IN=... OUT=..."
	@echo "  make bench-curves COST=<0/1> IN=<accepted.csv> OUT=<bench.csv>"
	@echo "  make bench-iso COST=<0/1> IN=<accepted.csv> OUT=<bench_iso.csv>"
	@echo ""
	@echo "Env for Sage: SMALL_BOUND=$(SMALL_BOUND), PRP=$(PRP), SKIP_TWIST=$(SKIP_TWIST), ALGO=$(ALGO)"
