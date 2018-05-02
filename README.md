# Sciluigi Workflows

Sciluigi workflows used to analyze data with multi-step execution

### Repository organization

The organization principle is as follows:
  * Every workflow is contained within its own executable Python script
  * Workflows load task classes via relative import

Therefore this repository contains a single folder which will contain all of
the task and workflow class definitions necessary. 

At some point in the future we may decide to wrap this up as an installable module,
but not yet.
