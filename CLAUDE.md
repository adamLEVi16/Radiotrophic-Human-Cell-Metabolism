# CLAUDE.md

## Workflow: Plan → Execute → Review

Don't default to running everything in one pass at max effort.

1. **Plan first.** Before touching files, state: the goal, files in/out of scope,
   constraints, edge cases, and what "done" looks like. Wait for approval on
   anything non-trivial or irreversible.
2. **Execute.** Do the mechanical work — edits, tests, boilerplate — without
   re-deriving the plan each step.
3. **Review.** Before declaring done, check the diff against the original plan.
   Flag anything overbuilt, fragile, or untested.

For genuinely large tasks (migrations, multi-file refactors), use subagents
with different models rather than one long session — e.g. a `plan` subagent
on a stronger model, an `execute` subagent on a faster/cheaper one. Configure
these in `.claude/agents/`, not as instructions here — CLAUDE.md can't switch
models mid-session on its own.

## Effort

- Default: medium — normal edits, standard debugging.
- High: architecture decisions, ambiguous planning, anything hard to undo.
- Don't max effort for formatting, renames, or boilerplate.

## Stop conditions

State one explicitly for any multi-step task before starting. Examples:
- Stop when tests pass and lint is clean.
- Stop after 3 failed attempts on the same approach; report the blocker instead
  of retrying with minor variations.
- Stop when the plan lists every file touched, every risk, and a rollback path.

Without a stop condition, don't start a loop — say so and ask for one.

## Constraints (edit per project)

- Keep tests passing.
- Preserve existing public interfaces unless explicitly told to break them.
- Match existing code style — don't reformat untouched code.
- Never touch secrets, credentials, or `.env` files.

## Handoff

When ending a session on unfinished work, write `memory/handoff.md`:

- Goal
- Completed
- Failed attempts (and why)
- Open blockers
- Files changed
- Next concrete action

Start the next session by reading that file first, not by re-reading the
whole prior conversation.

## Definition of done (for non-trivial tasks)

- [ ] Tests pass
- [ ] Lint passes
- [ ] Changed files summarized
- [ ] Risks called out
- [ ] Rollback path noted if the change is risky
