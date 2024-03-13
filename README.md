
# wdl-test-workflows
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)

This repository contains four basic workflows that can be used as introductory examples for users dipping their toes into the WDL waters, or equally as validation workflows for experienced users testing a new environment/configuration. 

# Usage

Test workflows are intended to be run in order of complexity:

- **helloHostname**: makes sure that your server is set up and the environment of the jobs is a valid working environment.
- **helloDockerHostname**: same as above but with the ubuntu:latest Docker container under the hood.
- **parseBatchFile**: tests access to a publicly available file in your filesystem, as well as the ability of the Cromwell server to parse that file and kick off a scatter of parallel jobs.
- **variantCalling**: tests whether the Cromwell server can do a multi-step, scientifically-relevant mini-workflow.

For Fred Hutch users that are new to WDL, we recommend using [PROOF](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/) to submit these workflows directly to the on-premise HPC cluster, as it simplifies interaction with Cromwell and provides a user-friendly front-end for job submission and tracking. For users outside of Fred Hutch or more advanced users who who would like to run the workflow locally, command line execution is relatively simple: 
```
java -jar cromwell-86.jar run helloHostname.wdl
java -jar cromwell-86.jar run helloDockerHostname.wdl
java -jar cromwell-86.jar run parseBatchFile.wdl --inputs parseBatchFile-inputs.json
java -jar cromwell-86.jar run variantCalling.wdl --inputs variantCalling-inputs.json
```
Although Cromwell is demonstrated here, these pipeline is not specific to Cromwell and can be run using whichever WDL execution method you prefer ([miniwdl](https://github.com/chanzuckerberg/miniwdl), [Terra](https://terra.bio/), [HealthOmics](https://docs.aws.amazon.com/omics/latest/dev/workflows.html), etc.).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on our [issue tracker](https://github.com/getwilds/wdl-test-workflows/issues).

## Contributing

If you would like to contribute to this WILDS WDL workflow, see our [contribution guidelines](.github/CONTRIBUTING.md) as well out our [WILDS Contributor Guide](https://getwilds.org/guide/) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.

