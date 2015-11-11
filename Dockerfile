# from the base image
FROM andrewosh/binder-base

# switching to main user
USER main

# adding requirements for my environment
ADD environment.yml environment.yml
