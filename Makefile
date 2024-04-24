# https://www.gnu.org/software/make/manual/html_node/Automatic-Variables.html#Automatic-Variables



HOT_MODULES =                                       \
	.                                               \

MODULES =                                           \
	$(HOT_MODULES)


OTHER_MODULES =                                     \

ALL_MODULES =                                       \
	$(OTHER_MODULES)                                \
	$(MODULES)


.PHONY: all

all: status

# Hot modules
status: $(addprefix STATUS_,$(HOT_MODULES))
add: $(addprefix ADD_,$(HOT_MODULES))
commit: $(addprefix COMMIT_,$(HOT_MODULES))
commit-no-edit: check_commit_msg $(addprefix COMMIT_,$(HOT_MODULES))
pull: $(addprefix PULL_,$(HOT_MODULES))
push: $(addprefix PUSH_,$(HOT_MODULES))

# All own modules
Status: $(addprefix STATUS_,$(MODULES))
Add: $(addprefix ADD_,$(MODULES))
Commit: check_commit_msg $(addprefix COMMIT_,$(MODULES))
Pull: $(addprefix PULL_,$(MODULES))
Push: $(addprefix PUSH_,$(MODULES))

# All modules: own and other external modules
# CLONE: $(addprefix CLONE_,$(ALL_MODULES))
PULL: $(ALL_MODULES) $(addprefix PULL_,$(ALL_MODULES))
REMOTE_GET-URL: $(addprefix REMOTE_GET-URL_,$(ALL_MODULES))

# Static pattern rules
$(addprefix ADD_,$(MODULES)): ADD_%: %
	@step "git -C $< add -A"
	try git -C $< add -A
	@next

# Static pattern rules
$(addprefix STATUS_,$(MODULES)): STATUS_%: %
	@echo ' '
	@echo "Status of $< " | sed -r 's/ \. / Root Repository /g'
	@echo "──────────────────────────────────────────────────────────────── "
	@git -C $< status

# Static pattern rules
$(addprefix COMMIT_,$(MODULES)): COMMIT_%: %
	@step "Trying to commit on $< with message \"$(m)\"" | sed -r 's/ \. / Root Repository /g'
	@echo "\n──────────────────────────────────────────────────────────────── "
	@echo ' '
	@try git -C $< commit -a $(m)
	@next

# Static pattern rules
$(addprefix PULL_,$(ALL_MODULES)): PULL_%: %
	@step "Trying to pull $< " | sed -r 's/ \. / Root Repository /g'
	@echo "\n──────────────────────────────────────────────────────────────── "
	@echo ' '
	@if git -C $< pull origin main ; then : ; else git -C $< push origin main ; fi
	@next

# Static pattern rules
$(addprefix PUSH_,$(MODULES)): PUSH_%: %
	@step "Trying to push $< " | sed -r 's/ \. / Root Repository /g'
	@echo "\n──────────────────────────────────────────────────────────────── "
	@echo ' '
	@if git -C $< push origin main ; then : ; else git -C $< push origin main ; fi
	@next


$(addprefix REMOTE_GET-URL_,$(ALL_MODULES)): REMOTE_GET-URL_%: %
	@echo -n "[\"$<\"]=\""
	@echo -n `git -C $< remote get-url origin`
	@echo -n "\"\n"


# 2020.12.12. (december 12, szombat), 14:54
check_commit_msg:
ifeq ($(m),)
	# $(warning No commit message defined)
	$(eval m := trivial)
	# $(info Automatic commit message is "$(m)")
	$(eval m := -m $(m))
else
	# $(info Your commit message is "$(m)")
	$(eval m := -m $(m))
endif


dummy:


show:
	git show --color --pretty=format:%b
