# ratesb: Rate Law Analysis for SBML and Antimony Models

Welcome to the web version of `ratesb`, a robust platform designed to analyze rate laws in SBML and Antimony models straight from your browser! This repository houses the web application's codebase, which offers the same core functionality as the `ratesb_python` package. (package will be released soon)

## Features

- Analyze rate laws in SBML and Antimony models via a user-friendly web interface.
- View a comprehensive list of warnings and errors in real-time.
- Add customized rate law and detect if your model is following the rate law forms.

## Standard Rate Law List

The standard rate law list encompasses a collection of commonly used rate laws, including:

- Zeroth order
- Uni-directional mass action
- Uni-term with the moderator
- Bi-directional mass action
- Bi-terms with the moderator
- Michaelis–Menten kinetics without explicit enzyme
- Michaelis–Menten kinetics with an explicit enzyme
- Hill equation
- Fraction format other than MM, MMCAT, and HILL

Each rate law is represented with specific attributes like the expression, optional symbols, and power-limited species, similar to the custom rate law classifications used in the `ratesbpython` package. These rate laws serve as a benchmark for the analysis and ensure a consistent user experience.

## Development Guide

### Prerequisites

Ensure you have the following software installed:

- Node.js
- npm
- Any other dependencies listed in `package.json`.

### Setup

1. Clone the repository:
```
git clone https://github.com/sys-bio/ratesb.git
```

2. Navigate to the repository directory:
```
cd ratesb
```

3. Install the necessary dependencies:
```
npm install
```

### Running the Application Locally

To run the application in a local development environment:

```
npx live-server
```

## Deployment Using GitHub Pages

`ratesb` is deployed using GitHub Pages. To deploy:

1. **Push your code to GitHub**: Ensure your latest codebase is pushed to your GitHub repository.

2. **Automatic Deployment Workflow**: Once the code is pushed, the GitHub Pages deployment workflow will automatically trigger and start the deployment process.

3. **Caution**: Avoid creating multiple deployment workflows in a short span as this might lead to potential crashes or overlaps in the workflow.

That's it! Your updates will reflect on the GitHub Pages-hosted site shortly after the deployment workflow completes.

## Contributing

We're always open to contributions! Whether it's bug reports, feature requests, or code tweaks, your inputs are invaluable. Please make a pull request or open an issue on our GitHub page.

## License

`ratesb` operates under the MIT license. Check the LICENSE file for further information.

## Future Improvements

We aim to continuously refine the performance of classifying Michaelis-Menten rate laws by evolving the underlying algorithm.

## Contact

For further inquiries or feedback, reach out to:

Longxuan Fan: longxuan@usc.edu