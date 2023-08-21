/** Provides Information on SBML Kinetics Laws */

class KineticLaw {
    constructor( libsbmlModel, libsbmlReaction, libsbmlKinetics, reaction, functionDefinitions, pyodide, processSBML, checks, customClassificationList) {
        this.libsbmlModel = libsbmlModel;
        this.libsbmlReaction = libsbmlReaction;
        this.pyodide = pyodide;
        this.libsbmlKinetics = libsbmlKinetics;
        this.formula = this.libsbmlKinetics.getFormula();
        this.reaction = reaction;
        this.model = processSBML;
        this.customClassificationList = customClassificationList;
        this._curDepth = 0;
        try {
            this.symbols = this._getSymbols();
        } catch (e) {
            this.symbols = [];
        }
        var symbolsLength = this.symbols.length;
        for (var i = 0; i < symbolsLength; i++) {
            if (typeof this.symbols[i] !== 'string' || this.symbols[i] === "") {
                this.symbols.splice(i, 1);
            }
            i--;
            symbolsLength--;
        }
        if (functionDefinitions == null) {
            this.expandedFormula = formula;
        } else {
            this.expandFormula(functionDefinitions);
        }
        this.sympyFormula = this.expandedFormula.replaceAll("^", "**");
        this.classificationCp = {};
        this.customClassification = [];
        this.analysis = {};
        this.checks = checks;
        if (this.formula.replace(/\s/g, '').length) {
            this.analysis["emptyRateLaw"] = false;
            this.classify();
        } else {
            this.analysis["emptyRateLaw"] = true;
        }
    }

    toString() {
        if (this.expandedFormula == null) {
            return self.formula;
        }
        return self.expandedFormula;
    }

    expandFormula(functionDefinitions) {
        var maxRecursion = 5;
        this.expandedFormula = this._expandFormula(this.formula, functionDefinitions, 0, maxRecursion);
    }

    classify() {
        var speciesInKineticLaw = [];
        var parametersInKineticLawOnly = [];
        var compartmentInKineticLaw = [];
        var othersInKineticLaw = [];
        var idsList = [...new Set(this.symbols)];
        var i;
        for (i = 0; i < idsList.length; i++) {
            if (this.model.speciesList.includes(idsList[i])) {
                speciesInKineticLaw.push(idsList[i]);
            } else if (this.model.parameterList.includes(idsList[i])) {
                parametersInKineticLawOnly.push(idsList[i]);
            } else if (this.model.compartmentList.includes(idsList[i])) {
                compartmentInKineticLaw.push(idsList[i]);
                othersInKineticLaw.push(idsList[i]);
            } else {
                othersInKineticLaw.push(idsList[i]);
            }
        }
        var parametersInKineticLaw = parametersInKineticLawOnly.concat(othersInKineticLaw);
        if (this.reaction.reactantList.length != 0) {
            idsList = idsList.concat(this.reaction.reactantList);
        }
        if (this.reaction.productList.length != 0) {
            idsList = idsList.concat(this.reaction.productList);
        }
        var idsList = [...new Set(idsList)];
        var symbolsSpace = idsList.join(" ");
        
        console.log(this.expandedFormula);
        if (idsList.length > 0) {
            this.pyodide.runPython(`
                ${idsList.toString()} = sympy.symbols("${symbolsSpace}")
            `);
        }
        this.pyodide.runPython(`
            try:
                kineticsSim = str(sympy.simplify(${this.expandedFormula}))
            except:
                kineticsSim = ${this.expandedFormula}
            js.kineticsSim = kineticsSim
        `);
        console.log(kineticsSim);
        if (isNaN(kineticsSim)) {
            this.analysis["pureNumber"] = false;
            var reactantList = this.reaction.reactantList;
            var productList = this.reaction.productList;
            var kinetics = this.expandedFormula;
    
            this.classificationCp["zerothOrder"] = this.isZerothOrder(speciesInKineticLaw);
            this.classificationCp["powerTerms"] = this.isPowerTerms(kinetics, kineticsSim);
            this.classificationCp["UNDR"] = this.isUNDR(reactantList, kinetics, kineticsSim, speciesInKineticLaw, idsList);
            this.classificationCp["UNMO"] = this.isUNMO(reactantList, kinetics, kineticsSim, speciesInKineticLaw, idsList);
            this.classificationCp["BIDR"] = this.isBIDR(reactantList, productList, kinetics, kineticsSim, speciesInKineticLaw, idsList);
            this.classificationCp["BIMO"] = this.isBIMO(reactantList, productList, kinetics, kineticsSim, speciesInKineticLaw, idsList);
            this.classificationCp["MM"] = this.isMM(kinetics, kineticsSim, idsList, speciesInKineticLaw, parametersInKineticLaw, reactantList);
            this.classificationCp["MMcat"] = this.isMMcat(kinetics, kineticsSim, idsList, speciesInKineticLaw, parametersInKineticLaw, reactantList, productList);
            this.classificationCp["Hill"] = this.isHill(kineticsSim, idsList, speciesInKineticLaw);
            this.classificationCp["Fraction"] = this.isFraction(kineticsSim, idsList, speciesInKineticLaw);
            this.classificationCp["Polynomial"] = this.isPolynomial(kineticsSim, idsList, speciesInKineticLaw);
    
            this.customClassification = this.customClassify(reactantList, productList, kineticsSim, idsList, speciesInKineticLaw, parametersInKineticLawOnly, compartmentInKineticLaw);
            console.log(this.classificationCp);
            console.log(this.customClassification);
            this.kineticLawAnalysis(speciesInKineticLaw, parametersInKineticLawOnly, compartmentInKineticLaw, othersInKineticLaw, idsList, kinetics, kineticsSim, reactantList, productList);
        } else {
            this.analysis["pureNumber"] = true;
        }
    }

    customClassify(reactantList, productList, kineticsSim, idsList, speciesInKineticLaw, parametersInKineticLawOnly, compartmentInKineticLaw) {
        function permute(arr) {
            // Base case: if the array has only one element, return it as a single-item permutation
            if (arr.length === 1) {
              return [arr];
            }
          
            // Array to store all permutations
            const permutations = [];
          
            // Iterate over each element in the array
            arr.forEach((element, index) => {
              // Create a copy of the array without the current element
              const remaining = [...arr.slice(0, index), ...arr.slice(index + 1)];
          
              // Generate permutations of the remaining elements
              const subPermutations = permute(remaining);
          
              // Add the current element to the beginning of each sub-permutation
              const mappedPermutations = subPermutations.map(subPermutation => [element, ...subPermutation]);
          
              // Add the mapped permutations to the main list
              permutations.push(...mappedPermutations);
            });
          
            return permutations;
        }

        function replaceOccurrences(reactantsInKineticLaw, productsInKineticLaw, enzymeList, compartmentInKineticLaw, parametersInKineticLawOnly, kineticsSim) {
            const permutedReactants = permute(reactantsInKineticLaw);
            const permutedProducts = permute(productsInKineticLaw);
            if (permutedReactants.length === 0) {
                permutedReactants.push([]);
            }
            if (permutedProducts.length === 0) {
                permutedProducts.push([]);
            }
            console.log(permutedReactants)
            console.log(permutedProducts)
            const ret = [];
            permutedReactants.forEach(reactantPerm => {
                permutedProducts.forEach(productPerm => {
                    const symbolPattern = /[^\w_]/g;
                    const symbols = kineticsSim.split(symbolPattern);
                    
                    const replacedSymbols = symbols.map(symbol => {
                        if (reactantPerm.includes(symbol)) {
                            var index = reactantPerm.indexOf(symbol) + 1;
                            return 'reactant' + index;
                        } else if (productPerm.includes(symbol)) {
                            var index = productPerm.indexOf(symbol) + 1;
                            return 'product' + index;
                        } else if (enzymeList.includes(symbol)) {
                            return 'enzyme';
                        } else if (compartmentInKineticLaw.includes(symbol)) {
                            return 'compartment';
                        } else if (parametersInKineticLawOnly.includes(symbol)) {
                            return 'parameter';
                        } else {
                            return symbol;
                        }
                    });
                
                    const nonAlphanumericChars = kineticsSim.match(symbolPattern) || [];
                    let replacedString = '';
                    
                    replacedSymbols.forEach((symbol, index) => {
                        replacedString += symbol + (nonAlphanumericChars[index] || '');
                    });
                    ret.push(replacedString);
                });
            });
            
            return ret;
        }
        const reactantsInKineticLaw = [];
        const productsInKineticLaw = [];
        console.log(speciesInKineticLaw)
        console.log(productList)
        speciesInKineticLaw.forEach(species => {
            if (reactantList.includes(species)) {
                reactantsInKineticLaw.push(species);
            } else if (productList.includes(species)) {
                productsInKineticLaw.push(species);
            }
        });
        const ret = [];
        const enzymeList = speciesInKineticLaw.filter(element => 
            !reactantList.includes(element) && !productList.includes(element)
        );
        const replacedKineticsList = replaceOccurrences(reactantsInKineticLaw, productsInKineticLaw, enzymeList, compartmentInKineticLaw, parametersInKineticLawOnly, kineticsSim);
        console.log(replacedKineticsList)
        this.customClassificationList.forEach(item => {
            const kineticsExpression = item.expression.replaceAll("^", "**");;
            const optionalSymbols = item.optional_symbols;
            const powerLimitedSpecies = item.power_limited_species;

            var classifiedTrue = false;
            for (let i = 0; i < replacedKineticsList.length; i++) {
                const replacedKinetics = replacedKineticsList[i];
                // Replace symbol names with their actual values in Python
                let pythonCode = `
def lower_powers(expr: sympy.Expr, keep=[]):
    """
    This function lowers the powers of a sympy expression to 1.
    It does not lower the powers of symbols specified in the "keep" list,
    and it only lowers positive integer powers.

    Parameters:
    expr (sympy.Expr): The sympy expression to transform.
    keep (list of sympy.Symbol): A list of symbols whose powers should not be lowered.

    Returns:
    sympy.Expr: The transformed sympy expression.
    """

    # Function to replace each instance of Pow in the expression
    def replace_if_applicable(base, exp):
        # If the exponent is a positive integer and the base is a symbol
        # not in the "keep" list, replace the Pow with the base
        if exp.is_integer and exp > 1 and base.is_symbol and base not in keep:
            return base
        else:
            # Otherwise, keep the Pow as it is
            return sympy.Pow(base, exp)

    # Use the replace method to apply the function to each Pow in the expression
    return expr.replace(sympy.Pow, replace_if_applicable)

def get_all_expr(expr, optional_symbols):
    """
    This function generates all possible expressions by optionally including
    each symbol in optional_symbols in the given expression.

    Parameters:
    expr (sympy.Expr): The sympy expression to transform.
    optional_symbols (list of sympy.Symbol): A list of symbols that can be optionally included.

    Returns:
    list of sympy.Expr: A list of all possible sympy expressions.
    """
    # Generate all possible combinations of the optional_symbols
    all_combinations = list(chain(*map(lambda x: combinations(optional_symbols, x), range(0, len(optional_symbols)+1))))

    all_expr = []
    for combo in all_combinations:
        # Copy the expression
        temp_expr = expr
        # Replace all non-included symbols with 1
        for sym in optional_symbols:
            if sym not in combo:
                temp_expr = temp_expr.subs(sym, 1)
        all_expr.append(sympy.simplify(temp_expr))

    return all_expr
reactant1, reactant2, reactant3, product1, product2, product3, enzyme, compartment, parameter = sympy.symbols('reactant1 reactant2 reactant3 product1 product2 product3 enzyme compartment parameter')
kineticsExpression = ${kineticsExpression}
all_expr = get_all_expr(kineticsExpression, [${optionalSymbols.toString()}])
js.allExpr = str(all_expr)
lowered_rate = lower_powers(${replacedKinetics}, [${powerLimitedSpecies.toString()}])
js.comparisonResult = any(expr == sympy.simplify(lowered_rate) for expr in all_expr)
                `;
                
                console.log(`comparing ${replacedKinetics} with ${kineticsExpression}`)
        
                // Run Python code in Pyodide
                this.pyodide.runPython(pythonCode);
                console.log(allExpr)
                if (comparisonResult) {
                    ret.push({
                        name: item.name,
                        comparisonResult: true
                    });
                    classifiedTrue = true;
                    break;
                }
            };
            if (!classifiedTrue) {
                ret.push({
                    name: item.name,
                    comparisonResult: false
                });
            }
        });
        return ret;
    }    

    kineticLawAnalysis(speciesInKineticLaw, parametersInKineticLawOnly, compartmentInKineticLaw, othersInKineticLaw, idsList, kinetics, kineticsSim, reactantList, productList) {
        this.assignPythonVals();
        const floatingSpecies = this.checkFloatingSpecies(speciesInKineticLaw, reactantList);
        this.analysis["floatingSpecies"] = floatingSpecies[0];
        this.analysis["boundaryFloatingSpecies"] = floatingSpecies[1];
        if (this.checks["reversibilityCheck"]) {
            this.analysis["inconsistentProducts"] = this.checkIrreversibleProduct(speciesInKineticLaw, productList);
        } else {
            this.analysis["inconsistentProducts"] = [];
        }
        this.analysis["nonIncreasingSpecies"] = this.checkNonIncreasingSpecies(speciesInKineticLaw, reactantList, productList, kinetics);
        if (this.checks["namingConvention"]) {
            this.analysis["namingConvention"] = this.checkNamingConventions(parametersInKineticLawOnly, kineticsSim, idsList);
        } else {
            this.analysis["namingConvention"] = {'k': [], 'K': [], 'V': []};
        }
        if (this.checks["formattingConvention"]) {
            this.analysis["formattingConvention"] = this.checkFormattingConventions(kinetics, kineticsSim, idsList, compartmentInKineticLaw, parametersInKineticLawOnly, speciesInKineticLaw);
        } else {
            this.analysis["formattingConvention"] = 0;
        }
        this.analysis["annotation"] = this.checkAnnotation();
        this.analysis["nonConstParams"] = this.constParamCheck(parametersInKineticLawOnly);
        this.assignPythonSymbols();
    }

    assignPythonVals() {
        this.analysis["uninitializedParameters"] = [];
        this.analysis["uninitializedCompartments"] = [];
        var object;
        var val;
        var id;
        for (var i = 0; i < this.libsbmlModel.getNumParameters(); i++) {
            object = this.libsbmlModel.getParameter(i);
            val = object.getValue();
            id = object.getId();
            if (val == 0 || isNaN(val)) {
                // need to give warning if val is not assigned
                val = 1;
                this.analysis["uninitializedParameters"].push(id);
            }
            this.pyodide.runPython(`
                ${id} = ${val}
            `);
        }
        for (var i = 0; i < this.libsbmlModel.getNumCompartments(); i++) {
            object = this.libsbmlModel.getCompartment(i);
            val = object.getSize();
            id = object.getId();
            if (val == 0 || isNaN(val)) {
                val = 1;
                this.analysis["uninitializedCompartments"].push(id);
            }
            this.pyodide.runPython(`
                ${id} = ${val}
            `);
        }
    }

    assignPythonSymbols() {
        var id;
        for (var i = 0; i < this.libsbmlModel.getNumParameters(); i++) {
            id = this.libsbmlModel.getParameter(i).getId();
            this.pyodide.runPython(`
                ${id} = sympy.symbols("${id}")
            `);
        }
        for (var i = 0; i < this.libsbmlModel.getNumCompartments(); i++) {
            id = this.libsbmlModel.getCompartment(i).getId();
            this.pyodide.runPython(`
                ${id} = sympy.symbols("${id}")
            `);
        }
    }

    /**
     * Check if the rate law contains all reactants
     * If reversible, check if the rate law contains all products
     * Return a list of species that violates this rule
     * @param {[]} speciesInKineticLaw list-species in the kinetics
     * @param {[]} reactantList list-reactants of the reaction
     * @return {[]} list of species violating the rule
     */
    checkFloatingSpecies(speciesInKineticLaw, reactantList) {
        var floatingSpecies = [];
        var boundaryFloatingSpecies = [];
        reactantList.forEach(reactant => {
            if (!speciesInKineticLaw.includes(reactant)) {
                if (this.reaction.boundarySpeciesList.includes(reactant)) {
                    boundaryFloatingSpecies.push(reactant);
                } else {
                    floatingSpecies.push(reactant);
                }
            }
        });
        return [floatingSpecies, boundaryFloatingSpecies];
    }

    checkIrreversibleProduct(speciesInKineticLaw, productList) {
        var inconsistentProducts = [];
        if (!this.libsbmlReaction.getReversible()) {
            productList.forEach(product => {
                if (speciesInKineticLaw.includes(product)) {
                    inconsistentProducts.push(product);
                }
            });
        }
        return inconsistentProducts;
    }

    /**
     * Check if the rate law is increasing as each reactant increases
     * If reversible, check if the rate law is decreasing as each product increases
     * Return a list of species that violates this rule
     * @param {[]} speciesInKineticLaw list-species in the kinetics
     * @param {[]} reactantList list-reactants of the reaction
     * @param {[]} productList list-products of the reaction
     * @param {string} kinetics string-kinetics
     * @return {[[], []]} list of species violating the rule, first list contains the reactants, second list contains the products.
     */
    checkNonIncreasingSpecies(speciesInKineticLaw, reactantList, productList, kinetics) {
        var nonIncreasingSpecies = [[],[]];

        this._assignPythonLocalVals();
        var prev;
        reactantList.forEach(reactant => {
            this._assignOtherSpeciesVals(reactant, speciesInKineticLaw, 1);
            prev = -Math.Infinity;
            for (var i = -0.9; i <= 0.9; i += 0.1) {
                this.pyodide.runPython(`
                    ${reactant} = ${Math.pow(10, i)}
                    js.curr = ${kinetics}
                `);
                if (curr <= prev) {
                    nonIncreasingSpecies[0].push(reactant);
                    break;
                }
                prev = curr;
            }
            this._assignPythonLocalSymbols(speciesInKineticLaw);
        });
        if (this.libsbmlReaction.getReversible()) {
            productList.forEach(product => {
                prev = Math.Infinity;
                if (speciesInKineticLaw.includes(product)) {
                    this._assignOtherSpeciesVals(product, speciesInKineticLaw, 1);
                    for (var i = -0.9; i <= 0.9; i += 0.1) {
                        this.pyodide.runPython(`
                            ${product} = ${Math.pow(10, i)}
                            js.curr = ${kinetics}
                        `);
                        if (curr >= prev) {
                            nonIncreasingSpecies[1].push(product);
                            break;
                        }
                        prev = curr;
                    }
                    this._assignSpeciesSymbols(speciesInKineticLaw);
                }
            });
        }
        this._assignPythonLocalSymbols();
        return nonIncreasingSpecies;
    }

    checkNamingConventions(parametersInKineticLaw, kineticsSim, idsList) {
        var ret = {'k': [], 'K': [], 'V': []};
        if (this.classificationCp["zerothOrder"]) {
            ret['k'] = this._checkSymbolsStartWith('k', parametersInKineticLaw);
        } else if (this.classificationCp["powerTerms"]) {

        } else if (this.classificationCp["UNDR"]) {
            ret['k'] = this._checkSymbolsStartWith('k', parametersInKineticLaw);
        } else if (this.classificationCp["UNMO"]) {
            ret['k'] = this._checkSymbolsStartWith('k', parametersInKineticLaw);
        } else if (this.classificationCp["BIDR"]) {
            ret['k'] = this._checkSymbolsStartWith('k', parametersInKineticLaw);
        } else if (this.classificationCp["BIMO"]) {
            ret['k'] = this._checkSymbolsStartWith('k', parametersInKineticLaw);
        } else if (this.classificationCp["MM"]) {
            var eq = this._numeratorDenominator(kineticsSim, idsList);
            var eq0 = [];
            var eq1 = [];
            parametersInKineticLaw.forEach(param => {
                if (eq[0].includes(param)) {
                    eq0.push(param);
                } else if (eq[1].includes(param)) {
                    eq1.push(param);
                }
            });
            ret['V'] = this._checkSymbolsStartWith('V', eq0);
            ret['K'] = this._checkSymbolsStartWith('K', eq1);
        } else if (this.classificationCp["MMcat"]) {
            var eq = this._numeratorDenominator(kineticsSim, idsList);
            var eq0 = [];
            var eq1 = [];
            parametersInKineticLaw.forEach(param => {
                if (eq[0].includes(param)) {
                    eq0.push(param);
                } else if (eq[1].includes(param)) {
                    eq1.push(param);
                }
            });
            ret['K'] = this._checkSymbolsStartWith('K', eq0);
            ret['K'] = this._checkSymbolsStartWith('K', eq1);
        } else if (this.classificationCp["Hill"]) {
            var eq = this._numeratorDenominator(kineticsSim, idsList);
            var eq0 = [];
            var eq1 = [];
            parametersInKineticLaw.forEach(param => {
                if (eq[0].includes(param)) {
                    eq0.push(param);
                } else if (eq[1].includes(param)) {
                    eq1.push(param);
                }
            });
            ret['K'] = this._checkSymbolsStartWith('K', eq1);
        } else if (this.classificationCp["Fraction"]) {
            
        } else if (this.classificationCp["Polynomial"]) {
            
        }
        return ret;
    }

    /**
     * Checks the formatting conventions of a rate law expression based on compartments,
     * parameters, and species. Handles specific cases for MM, MMcat, Hill, and Fraction classifications.
     * 
     * @param {string} kinetics - The original rate law to analyze.
     * @param {string} kineticsSim - The simplified rate law to analyze.
     * @param {string[]} idsList - An array of element IDs to consider.
     * @param {string[]} compartmentInKineticLaw - An array of compartments in the rate law.
     * @param {string[]} parametersInKineticLawOnly - An array of parameters in the rate law.
     * @param {string[]} speciesInKineticLaw - An array of species in the rate law.
     * @returns {number} A number representing the formatting convention status, where:
     *                   0 = formatted correctly,
     *                   1 = elements not in alphabetical order,
     *                   2 = formatting convention not followed,
     *                   >= 3 = rate law in fraction form,
     *                   3 = denominator not in alphabetical order,
     *                   4 = numerator and denominator not in alphabetical order,
     *                   5 = numerator convention not followed and denominator not in alphabetical order,
     *                   6 = denominator convention not followed,
     *                   7 = numerator not in alphabetical order and denominator convention not followed,
     *                   8 = numerator and denominator convention not followed.
     */
    checkFormattingConventions(kinetics, kineticsSim, idsList, compartmentInKineticLaw, parametersInKineticLawOnly, speciesInKineticLaw) {
        var ret = 0;
        if (this.classificationCp["MM"] ||
            this.classificationCp["MMcat"] ||
            this.classificationCp["Hill"] ||
            this.classificationCp["Fraction"]
            ) {
            var eq = this.splitFraction(kinetics);
            var eq0 = eq[0];
            var eq1 = eq[1];
            ret = this._checkExpressionFormat(eq0, compartmentInKineticLaw, parametersInKineticLawOnly);
            ret += 3 * this._checkExpressionFormat(eq1, compartmentInKineticLaw, parametersInKineticLawOnly);
        } else {
            ret = this._checkExpressionFormat(kinetics, compartmentInKineticLaw, parametersInKineticLawOnly);
        }
        return ret;
    }

    splitFraction(fractionString) {
        // Split the fraction at the '/' character
        let splitFraction = fractionString.split('/');
    
        // Check if there's a numerator and a denominator
        if (splitFraction.length !== 2) {
            return ['',''];
        }
    
        // Return the numerator and the denominator as an object
        return [splitFraction[0].trim(), splitFraction[1].trim()];
    }

    /**
     * This function checks the SBOTerm annotation against various sets of 
     * acceptable terms for different rate law types in systems biology.
     *
     * The function checks whether the term belongs to a particular set of SBO terms, 
     * representing different rate law types (UNDR, UNMO, BIDR, BIMO, MM, MMCAT). 
     * If the term doesn't belong to the corresponding set for the detected type, 
     * the function will return a unique warning code associated with each rate law type.
     * 
     * If the term is found within the respective set of terms, the function will return 0,
     * indicating no issues found.
     * 
     * @returns {number} - Returns 0 if the SBOTerm is in the corresponding set. 
     * Otherwise, returns a unique number representing the type of issue (1-5).
     */
    checkAnnotation() {
        const sboTerm = this.libsbmlKinetics.getSBOTerm();
        if (sboTerm < 0) {
            return 0;
        }
        if (this.classificationCp["UNDR"]) {
            const undrSBOs = [41, 43, 44, 45, 47, 49, 50, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 140, 141, 142, 143, 144, 145, 146, 163, 166, 333, 560, 561, 562, 563, 564, 430, 270, 458, 275, 273, 379, 440, 443, 451, 454, 456, 260, 271, 378, 387, 262, 265, 276, 441, 267, 274, 444, 452, 453, 455, 457, 386, 388, 442, 277, 445, 446, 447, 448, 266, 449, 450];
            if (!undrSBOs.includes(sboTerm)) {
                return 1;
            }
        } else if (this.classificationCp["UNMO"]) {
            const unmoSBOs = [41, 43, 44, 45, 47, 49, 50, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 140, 141, 142, 143, 144, 145, 146, 163, 166, 333, 560, 561, 562, 563, 564];
            if (!unmoSBOs.includes(sboTerm)) {
                return 2;
            }
        } else if (this.classificationCp["BIDR"] || this.classificationCp["BIMO"]) {
            const biSBOs = [42, 69, 78, 88, 109, 646, 70, 71, 74, 79, 80, 81, 84, 89, 99, 110, 120, 130, 72, 73, 75, 76, 77, 82, 83, 85, 86, 87, 90, 91, 92, 95, 100, 101, 102, 105, 111, 112, 113, 116, 121, 122, 123, 126, 131, 132, 133, 136, 93, 94, 96, 97, 98, 103, 104, 106, 107, 108, 114, 115, 117, 118, 119, 124, 125, 127, 128, 129, 134, 135, 137, 138, 139];
            if (!biSBOs.includes(sboTerm)) {
                return 3;
            }
        } else if (this.classificationCp["MM"]) {
            const mmSBOs = [28, 29, 30, 31, 199];
            if (!mmSBOs.includes(sboTerm)) {
                return 4;
            }
        } else if (this.classificationCp["MMcat"]) {
            const mmcatSBOs = [28, 29, 30, 31, 199, 430, 270, 458, 275, 273, 379, 440, 443, 451, 454, 456, 260, 271, 378, 387, 262, 265, 276, 441, 267, 274, 444, 452, 453, 455, 457, 386, 388, 442, 277, 445, 446, 447, 448, 266, 449, 450];
            if (!mmcatSBOs.includes(sboTerm)) {
                return 5;
            }
        } else if (this.classificationCp["Hill"]) {
            const mmcatSBOs = [192, 195, 198];
            if (!mmcatSBOs.includes(sboTerm)) {
                return 6;
            }
        }
        return 0;
    }

    constParamCheck(parametersInKineticLawOnly) {
        const nonConstParams = [];
        parametersInKineticLawOnly.forEach(param => {
            const libsbmlParam = this.libsbmlModel.getParameter(param);
            if (!libsbmlParam.getConstant()) {
                nonConstParams.push(param);
            }
        });
        return nonConstParams;
    }

    escapeRegex(string) {
        return string.replace(/[/\-\\^$*+?.()|[\]{}]/g, '\\$&');
    }

    /**
     * Finds the positions of elements in a given rate law and returns an object containing
     * the largest position, smallest position, and a flag indicating if the positions are increasing.
     * 
     * @param {string[]} elementList - An array of elements to find in the rate law.
     * @param {string} rateLaw - The rate law in which to search for the elements.
     * @returns {Object} An object containing the largest position, smallest position, and increasing flag.
     */
    _findPositionsInRateLaw(elementList, rateLaw) {
        var largestPosition = -1;
        var smallestPosition = Infinity;
        var prevPosition = -1;
        var increasingFlag = true;

        elementList.forEach(element => {
            var pattern = new RegExp("\\b" + this.escapeRegex(element) + "\\b");
            var index = rateLaw.search(pattern);
            if (index >= 0) {
                largestPosition = Math.max(index, largestPosition);
                smallestPosition = Math.min(index, smallestPosition);
                if (index < prevPosition) {
                    increasingFlag = false;
                }
                prevPosition = index;
            }
        });

        return {largestPosition, smallestPosition, increasingFlag};
    }

    /**
     * Checks the format of a rate law based on the positions of compartments, parameters, and species.
     * 
     * @param {string} kinetics - The rate law to analyze.
     * @param {string[]} compartmentInKineticLaw - An array of compartments in the rate law.
     * @param {string[]} parametersInKineticLawOnly - An array of parameters in the rate law.
     * @param {string[]} speciesInKineticLaw - An array of species in the rate law.
     * @returns {number} Returns 0 if formatted correctly, 1 if elements are not in alphabetical order, or 2 if the formatting convention is not followed.
     */
    _checkProductOfTermsFormat(kinetics, compartmentInKineticLaw, parametersInKineticLawOnly, speciesInKineticLaw) {
        const compStats = this._findPositionsInRateLaw(compartmentInKineticLaw, kinetics);
        const paramStats = this._findPositionsInRateLaw(parametersInKineticLawOnly, kinetics);
        const specStats = this._findPositionsInRateLaw(speciesInKineticLaw, kinetics);
        var isFormatted = (compStats.largestPosition < paramStats.smallestPosition || paramStats.smallestPosition === Infinity)
              && (paramStats.largestPosition < specStats.smallestPosition || specStats.smallestPosition === Infinity)
              && (compStats.largestPosition < specStats.smallestPosition || specStats.smallestPosition === Infinity);
        var increasingFlag = compStats.increasingFlag && paramStats.increasingFlag && specStats.increasingFlag;
        if (isFormatted) {
            if (increasingFlag) {
                return 0;
            } else {
                return 1;
            }
        }
        return 2;
    }


    /**
     * Analyzes a rate law expression and checks its format based on the positions of compartments,
     * parameters, and species.
     * 
     * @param {string} kinetics - The rate law to analyze.
     * @param {string[]} compartmentInKineticLaw - An array of compartments in the rate law.
     * @param {string[]} parametersInKineticLawOnly - An array of parameters in the rate law.
     * @returns {number} Returns 0 if formatted correctly, 1 if elements are not in alphabetical order, or 2 if the formatting convention is not followed.
     */
    _checkExpressionFormat(kinetics, compartmentInKineticLaw, parametersInKineticLawOnly) {
        compartmentInKineticLaw.sort();
        parametersInKineticLawOnly.sort();
        var productOfTerms = kinetics.split(/[+-]/);
        var ret = 0;
        productOfTerms.forEach(term => {
            var format = this._checkProductOfTermsFormat(term, compartmentInKineticLaw, parametersInKineticLawOnly, this.reaction.sortedSpeciesList);
            if (format > ret) {
                ret = format;
            }
        });
        return ret;
    }

    /**
     * check if symbols in the list start with given start
     * Case Sensitive
     * @param {String} start 
     * @param {[]} symbols 
     * @returns {[]} list of symbols not starting with given start
     */
    _checkSymbolsStartWith(start, symbols) {
        var ret = [];
        symbols.forEach(symbol => {
            if (!symbol.startsWith(start)) {
                ret.push(symbol);
            }
        });
        return ret;
    }

    _assignPythonLocalVals() {
        for (var i = 0; i < this.libsbmlKinetics.getNumParameters(); i++) {
            var param = this.libsbmlKinetics.getParameter(i);
            this.pyodide.runPython(`
                ${param.getId()} = ${param.getValue()}
            `);
        }
        for (var i = 0; i < this.libsbmlKinetics.getNumLocalParameters(); i++) {
            var param = this.libsbmlKinetics.getLocalParameter(i);
            this.pyodide.runPython(`
                ${param.getId()} = ${param.getValue()}
            `);
        }
    }

    _assignPythonLocalSymbols() {
        for (var i = 0; i < this.libsbmlKinetics.getNumParameters(); i++) {
            var param = this.libsbmlKinetics.getParameter(i).getId();
            this.pyodide.runPython(`
                ${param} = sympy.symbols("${param}")
            `);
        }
        for (var i = 0; i < this.libsbmlKinetics.getNumLocalParameters(); i++) {
            var param = this.libsbmlKinetics.getLocalParameter(i);
            this.pyodide.runPython(`
                ${param} = sympy.symbols("${param}")
            `);
        }
    }

    _assignOtherSpeciesVals(species, speciesList, val) {
        speciesList.forEach(otherSpecies => {
            if (otherSpecies !== species) {
                this.pyodide.runPython(`
                    ${otherSpecies} = ${val}
                `);
            }
        });
    }

    _assignSpeciesSymbols(speciesList) {
        speciesList.forEach(species => {
            this.pyodide.runPython(`
                ${species} = sympy.symbols("${species}")
            `);
        });
    }

    /**
     * Check whether the reaction with a kinetic law belongs to the type of Zeroth Order
     * Zeroth order classification rule: if there are no species in the kinetics
     * @param {[]} speciesInKineticLaw list-species in the kinetics
     * @returns {boolean} whether the rate law is zeroth order
     */
    isZerothOrder(speciesInKineticLaw) {
        return this._numSpeciesInKinetics(speciesInKineticLaw) == 0;
    }

    /**
     * Check whether the reaction with a kinetic law belongs to the type of Kinetics with power terms
     * Kinetics with power terms classification rule: if there is pow() or ** inside the kinetics,
     * except the pow(,-1) case as the possible Michaelis–Menten kinetics
     * @param {string} kinetics string-kinetics
     * @param {string} kineticsSim string-simplified kinetics
     * @returns {boolean} whether the rate law is the type of Kinetics with power terms
     */
    isPowerTerms(kinetics, kineticsSim) {
        return this._powerInKinetics(kinetics, kineticsSim);
    }

    /**
     * Tests whether the reaction belongs to the type of No Products
     * No products classification rule: if there are no products
     * @param {[]} productList list-products of the reaction
     * @returns {boolean} whether the reaction belongs to the type of No Products
     */
    isNoPrds(productList) {
        return this._numOfPrds(productList) == 0;
    }

    /**
     * Tests whether the reaction belongs to the type of Single Products
     * No products classification rule: if there is a single product
     * @param {[]} productList list-products of the reaction
     * @returns {boolean} whether the reaction belongs to the type of Single Products
     */
    isSinglePrd(productList) {
        return this._numOfPrds(productList) == 1;
    }

    /**
     * Tests whether the reaction belongs to the type of Double Products
     * No products classification rule: if there are two products
     * @param {[]} productList list-products of the reaction
     * @returns {boolean} whether the reaction belongs to the type of Double Products
     */
    isDoublePrds(productList) {
        return this._numOfPrds(productList) == 2;
    }

    /**
     * Tests whether the reaction belongs to the type of Multiple Products
     * No products classification rule: if there are more than two products
     * @param {[]} productList list-products of the reaction
     * @returns {boolean} whether the reaction belongs to the type of Multiple Products
     */
    isMulPrds(productList) {
        return this._numOfPrds(productList) > 2;
    }

    /**
     * Tests whether the reaction belongs to the type of No Reactants
     * No reactants classification rule: if there are no reactants
     * @param {[]} reactantList list-reactants of the reaction
     * @returns {boolean} whether the reaction belongs to the type of No Reactants
     */
    isNoRcts(reactantList) {
        return this._numOfRcts(reactantList) == 0;
    }

    /**
     * Tests whether the reaction belongs to the type of Single Reactants
     * No reactants classification rule: if there is a single product
     * @param {[]} reactantList list-reactants of the reaction
     * @returns {boolean} whether the reaction belongs to the type of Single Reactants
     */
    isSingleRct(reactantList) {
        return this._numOfRcts(reactantList) == 1;
    }

    /**
     * Tests whether the reaction belongs to the type of Double Reactants
     * No reactants classification rule: if there are two reactants
     * @param {[]} reactantList list-reactants of the reaction
     * @returns {boolean} whether the reaction belongs to the type of Double Reactants
     */
    isDoubleRcts(reactantList) {
        return this._numOfRcts(reactantList) == 2;
    }

    /**
     * Tests whether the reaction belongs to the type of Multiple Reactants
     * No reactants classification rule: if there are more than two reactants
     * @param {[]} reactantList list-reactants of the reaction
     * @returns {boolean} whether the reaction belongs to the type of Multiple Reactants
     */
    isMulRcts(reactantList) {
        return this._numOfRcts(reactantList) > 2;
    }

    /**
     * Tests whether the reaction belongs to the type of uni-directional mass reaction
     * Uni-directional mass reaction classification rule: 
     * 1) Kinetics is a single product of terms or with an additional const term
     * 2) The species inside the kinetics are only reactants
     * @param {[]} reactantList list-reactants of the reaction
     * @param {string} kinetics string-kinetics
     * @param {string} kineticsSim string-simplified kinetics
     * @param {[]} speciesInKineticLaw list-species in the kinetics
     * @param {[]} idsList list-id list including all the ids in kinetics, reactants and products
     * @returns {boolean} whether the rate law is UNDR
     */
    isUNDR(reactantList, kinetics, kineticsSim, speciesInKineticLaw, idsList) {
        var flag = false;
        if (this._isSingleProductOfTerms(kinetics, kineticsSim) && this._specsInKineticsAllRcts(speciesInKineticLaw, reactantList)) {
            flag = true;
        }
        if (speciesInKineticLaw.length == 1 && reactantList.length == 1) {
            if (speciesInKineticLaw[0] === reactantList[0]) {
                if (kinetics.split(speciesInKineticLaw[0]).length - 1 == 1) {
                    flag = true;
                } else if (kineticsSim.split(speciesInKineticLaw[0]).length - 1 == 1) {
                    flag = true;
                }
            }
        }
        try {
            var eq = this._numeratorDenominator(kineticsSim, idsList);
            for (var i = 0; i < speciesInKineticLaw.length; i++) {
                if (eq[1].includes(speciesInKineticLaw[i])) {
                    flag = false;
                }
            }
        } catch {

        }
        return flag;
    }

    /**
     * Tests whether the reaction belongs to the type of uni-term with moderator
     * Uni-term with moderator classification rule: 
     * 1) Kinetics is a single product of terms
     * 2) The species inside the kinetics are not only reactants
     * @param {[]} reactantList list-reactants of the reaction
     * @param {string} kinetics string-kinetics
     * @param {string} kineticsSim string-simplified kinetics
     * @param {[]} speciesInKineticLaw list-species in the kinetics
     * @param {[]} idsList list-id list including all the ids in kinetics, reactants and products
     * @returns {boolean} whether the rate law is UNMO
     */
    isUNMO(reactantList, kinetics, kineticsSim, speciesInKineticLaw, idsList) {
        var flag = false;
        if (this._isSingleProductOfTerms(kinetics, kineticsSim)
            && !this._specsInKineticsAllRcts(speciesInKineticLaw, reactantList)
            && this._numSpeciesInKinetics(speciesInKineticLaw) != 0) {
            flag = true;
        }
        if (speciesInKineticLaw.length == 1 && !compareLists(speciesInKineticLaw, reactantList) && !this._isDiffOfTwoProductsOfTerms(kinetics, kineticsSim)) {
            if (kinetics.split(speciesInKineticLaw[0]).length - 1 == 1) {
                flag = true;
            } else if (kineticsSim.split(speciesInKineticLaw[0]).length - 1 == 1) {
                flag = true;
            }
        }
        try {
            var eq = this._numeratorDenominator(kineticsSim, idsList);
            for (var i = 0; i < speciesInKineticLaw.length; i++) {
                if (eq[1].includes(speciesInKineticLaw[i])) {
                    flag = false;
                }
            }
        } catch {

        }
        return flag;
    }

    /**
     * Tests whether the reaction belongs to the type of bi-directional mass reaction
     * Bi-directional mass reactionclassification rule: 
     * 1) Kinetics is the difference of two product of terms, with species in each term
     * 2) The first term before - includes all the reactants
     *    while the second term after - includes all the products
     * @param {[]} reactantList list-reactants of the reaction
     * @param {[]} productList list-products of the reaction
     * @param {string} kinetics string-kinetics
     * @param {string} kineticsSim string-simplified kinetics
     * @param {[]} speciesInKineticLaw list-species in the kinetics
     * @param {[]} idsList list-id list including all the ids in kinetics, reactants and products
     * @returns {boolean} whether the rate law is BIDR
     */
    isBIDR(reactantList, productList, kinetics, kineticsSim, speciesInKineticLaw, idsList) {
        var flag = true;
        if (!this._isDiffOfTwoProductsOfTerms(kinetics, kineticsSim)) {
            flag = false;
        }
        if (!this._ProductOfTermsWithAllRctsOrPrds(kinetics, kineticsSim, speciesInKineticLaw, reactantList, productList)) {
            flag = false;
        }
        try {
            var eq = this._numeratorDenominator(kineticsSim, idsList);
            for (var i = 0; i < speciesInKineticLaw.length; i++) {
                if (eq[1].includes(speciesInKineticLaw[i])) {
                    flag = false;
                }
            }
        } catch {

        }
        return flag;
    }

    /**
     * Tests whether the reaction belongs to the type of bi-terms with moderator
     * Bi-terms with moderator classification rule: 
     * 1) Kinetics is the difference of two product of terms
     * 2) The first term before - does not include all the reactants
     *    while the second term after - does not include all the products
     * @param {[]} reactantList list-reactants of the reaction
     * @param {[]} productList list-products of the reaction
     * @param {string} kinetics string-kinetics
     * @param {string} kineticsSim string-simplified kinetics
     * @param {[]} speciesInKineticLaw list-species in the kinetics
     * @param {[]} idsList list-id list including all the ids in kinetics, reactants and products
     * @returns {boolean} whether the rate law is BIMO
     */
    isBIMO(reactantList, productList, kinetics, kineticsSim, speciesInKineticLaw, idsList) {
        var flag = true;
        if (this._numSpeciesInKinetics(speciesInKineticLaw) == 0) {
            flag = false;
        }
        if (!this._isDiffOfTwoProductsOfTerms(kinetics, kineticsSim)) {
            flag = false;
        }
        if (this._ProductOfTermsWithAllRctsOrPrds(kinetics, kineticsSim, speciesInKineticLaw, reactantList, productList)) {
            flag = false;
        }
        if (speciesInKineticLaw.length == 1 && reactantList.length == 1) {
            if (speciesInKineticLaw[0] === reactantList[0]) {
                if (kinetics.split(speciesInKineticLaw[0]).length - 1 == 1) {
                    flag = false;
                } else if (kineticsSim.split(speciesInKineticLaw[0]).length - 1 == 1) {
                    flag = false;
                }
            }
        }
        try {
            eq = this._numeratorDenominator(kineticsSim, idsList);
            if (speciesInKineticLaw.length > 0) {
                speciesInKineticLaw.forEach(species => {
                    if (eq[1].includes(species)) {
                        flag = false;
                    }
                })
            }
        } catch {

        }
        return flag;
    }

    isMM(kinetics, kineticsSim, idsList, speciesInKineticLaw, parametersInKineticLaw, reactantList) {
        var eq = this._numeratorDenominator(kineticsSim, idsList);
        var flagFr = false;
        speciesInKineticLaw.forEach(species => {
            if (eq[1].includes(species)) {
                flagFr = true;
            }
        });
        if (flagFr) {
            if (this._numSpeciesInKinetics(speciesInKineticLaw) == 1 && this._numOfRcts(reactantList) == 1) {
                if (this._MMSingleSpecInNumerator(kinetics, idsList, parametersInKineticLaw, reactantList)) {
                    return true;
                }
            }
        }
        return false;
    }

    isMMcat(kinetics, kineticsSim, idsList, speciesInKineticLaw, parametersInKineticLaw, reactantList, productList) {
        var eq = this._numeratorDenominator(kineticsSim, idsList);
        var flagFr = false;
        speciesInKineticLaw.forEach(species => {
            if (eq[1].includes(species)) {
                flagFr = true;
            }
        });
        if (flagFr) {
            if (this._numSpeciesInKinetics(speciesInKineticLaw) == 2 && this._numOfRcts(reactantList) == 1
                && !speciesInKineticLaw.some(item => productList.includes(item))) {
                if (this._MMSingleSpecInNumerator(kinetics, idsList, parametersInKineticLaw, reactantList)) {
                    return true;
                }
            }
        }
        return false;
    }

    isHill(kineticsSim, idsList, speciesInKineticLaw) {
        var eq = this._numeratorDenominator(kineticsSim, idsList);
        var flagFr = false;
        speciesInKineticLaw.forEach(species => {
            if (eq[1].includes(species)) {
                flagFr = true;
            }
        });
        if (flagFr) {
            if (this._numSpeciesInKinetics(speciesInKineticLaw) == 1) {
                if (this._HillFormat(kineticsSim, idsList, speciesInKineticLaw)) {
                    return true;
                }
            }
        }
        return false;
    }

    isFraction(kineticsSim, idsList, speciesInKineticLaw) {
        var flag = false;
        var eq = this._numeratorDenominator(kineticsSim, idsList);
        speciesInKineticLaw.forEach(species => {
            if (eq[1].includes(species)) {
                flag = true;
            }
        });
        return flag;
    }

    isPolynomial(kineticsSim, idsList, speciesInKineticLaw) {
        var flag = false;
        if (this._isPolynomial(kineticsSim, idsList) && speciesInKineticLaw.length > 0) {
            speciesInKineticLaw.forEach(species => {
                if (kineticsSim.includes(species)) {
                    flag = true;
                }
            });
        }
        return flag;
    }

    /**
     * Returns the number of species in the kinetic law
     * @param {[]} speciesInKineticLaw list-species in the kinetics
     * @return {int}
     */
    _numSpeciesInKinetics(speciesInKineticLaw) {
        return speciesInKineticLaw.length;
    }

    /**
     * Tests whether there is power term in the kinetic law: **, "pow", but not "pow(,-1)"
     * @param {string} kinetics string-kinetics
     * @param {string} kineticsSim string-simplified kinetics
     * @return {boolean}
     */
    _powerInKinetics(kinetics, kineticsSim) {
        if (kinetics.includes("pow(") && !kinetics.includes("-1)")) {
            return true;
        } else if (kinetics.includes("**")) {
            return true;
        } else if (kineticsSim.includes("**")) {
            return true;
        }
        return false;
    }

    /**
     * Returns the number of products in the reaction
     * @param {[]} productList list-products of the reaction
     * @return {int}
     */
    _numOfPrds(productList) {
        return productList.length;
    }

    /**
     * Returns the number of reactants in the reaction
     * @param {[]} reactantList list-reactants of the reaction
     * @return {int}
     */
    _numOfRcts(reactantList) {
        return reactantList.length;
    }

    /**
     * Tests whether the kinetics is a single product of terms
     * @param {string} kinetics string-kinetics
     * @param {string} kineticsSim string-simplified kinetics
     * @return {boolean}
     */
    _isSingleProductOfTerms(kinetics, kineticsSim) {
        var flag = true;
        if (kinetics.includes("+") || kinetics.includes("-")) {
            flag = false;
            if (kinetics.includes("e-") || kinetics.includes("exp(-")) {
                flag = true;
            }
        } else if (kineticsSim.includes("+") || kineticsSim.includes("-")) {
            flag = false;
            if (kinetics.includes("e-") || kineticsSim.includes("exp(-")) {
                flag = true;
            }
        }
        return flag;
    }

    /**
     * Tests whether all species in kinetics are reactants
     * @param {[]} speciesInKineticLaw list-species in the kinetics
     * @param {[]} reactantList list-reactants of the reaction
     * @return {boolean}
     */
    _specsInKineticsAllRcts(speciesInKineticLaw, reactantList) {
        if (reactantList.length > 0 && compareMaps(Counter(speciesInKineticLaw), Counter(reactantList))) {
            return true;
        }
        return false;
    }

    /**
     * Tests whether the kinetics is the difference between two product of terms
     * @param {string} kinetics string-kinetics
     * @param {string} kineticsSim string-simplified kinetics
     * @return {boolean}
     */
    _isDiffOfTwoProductsOfTerms(kinetics, kineticsSim) {
        var flag = false;
        if (!this._isSingleProductOfTerms(kinetics, kineticsSim)) {
            var terms = kinetics.split("-");
            if (terms.length == 2) {
                flag = true;
            }
        }
        return flag;
    }

    /**
     * Tests whether the kinetics with one/the other product terms with all reactants/products
     * @param {string} kinetics string-kinetics
     * @param {string} kineticsSim string-simplified kinetics
     * @param {[]} speciesInKineticLaw list-species in the kinetics
     * @param {[]} reactantList list-reactants of the reaction
     * @param {[]} productList list-products of the reaction
     * @return {boolean}
     */
    _ProductOfTermsWithAllRctsOrPrds(kinetics, kineticsSim, speciesInKineticLaw, reactantList, productList) {
        var flagKinetics = true;
        var flagKineticsSim = true;
        var terms = kinetics.split("-");
        if (terms.length == 2 && reactantList.length > 0 && productList.length > 0) {
            var term1 = terms[0];
            var term2 = terms[1];
            if (compareMaps(Counter(speciesInKineticLaw), Counter(reactantList.concat(productList)))) {
                reactantList.forEach(element => {
                    if (!term1.includes(element)) {
                        flagKinetics = false;
                    }
                });
                productList.forEach(element => {
                    if (!term2.includes(element)) {
                        flagKinetics = false;
                    }
                });
            } else {
                flagKinetics = false;
            }
        }
        terms = kineticsSim.split("-");
        if (terms.length == 2 && reactantList.length > 0 && productList.length > 0) {
            var term1 = terms[0];
            var term2 = terms[1];
            if (compareMaps(Counter(speciesInKineticLaw), Counter(reactantList.concat(productList)))) {
                reactantList.forEach(element => {
                    if (!term1.includes(element)) {
                        flagKineticsSim = false;
                    }
                });
                productList.forEach(element => {
                    if (!term2.includes(element)) {
                        flagKineticsSim = false;
                    }
                });
            } else {
                flagKineticsSim = false;
            }
        }
        return flagKinetics || flagKineticsSim;
    }

    /**
     * Tests whether kinetics is in the MM functional form with a single species in the numerator
     * @param {string} kinetics string-kinetics
     * @param {[]} idsList list-id list including all the ids in kinetics, reactants and products
     * @param {[]} parametersInKineticLaw list-parameters in the kinetics
     * @param {[]} reactantList list-reactants of the reaction
     * @return {boolean}
     */
    _MMSingleSpecInNumerator(kinetics, idsList, parametersInKineticLaw, reactantList) {
        var preSymbols = "";
        var i;
        for (i = 0; i < idsList.length; i++) {
            preSymbols += idsList[i];
            preSymbols += " ";
        }
        preSymbols = preSymbols.substring(0, preSymbols.length - 1);
        var preSymbolsComma = preSymbols.replaceAll(" ", ",");
        this.pyodide.runPython(`
                try:
                    ${preSymbolsComma} = sympy.symbols("${preSymbols}")
                    js.strangeFunc1 = 0
                except:
                    js.strangeFunc1 = 1
            `);
        this.pyodide.runPython(`
            try:
                expr = ${kinetics}
                js.strangeFunc2 = 0
            except:
                js.strangeFunc2 = 1
        `);
        if (strangeFunc1 + strangeFunc2 == 0) {
            if (parametersInKineticLaw.length >= 2) {
                for (var j = 0; j < parametersInKineticLaw.length; j++) {
                    for (var k = 0; k < parametersInKineticLaw.length; k++) {
                        for (var l = 0; l < parametersInKineticLaw.length; l++) {
                            for (var m = 0; m < parametersInKineticLaw.length; m++) {
                                if (k != j) {
                                    var preN = `${reactantList[0]} * ${parametersInKineticLaw[j]}`;
                                    var preD = `( ${reactantList[0]} + ${parametersInKineticLaw[k]} )`;
                                    var expr1Stat = `expr1 = ${preN}  / ${preD}`
                                    this.pyodide.runPython(`
                                        ${expr1Stat}
                                        if sympy.simplify(expr1) == sympy.simplify(expr):
                                            js.flag = 1
                                    `);
                                    if (typeof(flag) !== undefined) {
                                        return true;
                                    }
                                    if (parametersInKineticLaw.length >= 3) {
                                        if (l != j && l != k) {
                                            preN += ` * ${parametersInKineticLaw[l]}`
                                            expr1Stat = `expr1 = ${preN} / ${preD}`;
                                            this.pyodide.runPython(`
                                                ${expr1Stat}
                                                if sympy.simplify(expr1) == sympy.simplify(expr):
                                                    js.flag = 1
                                            `);
                                            if (typeof(flag) !== undefined) {
                                                return true;
                                            }
                                        }
                                        if (parametersInKineticLaw.length >= 4) {
                                            if (m != j && m != k && m != l) {
                                                preN += ` * ${parametersInKineticLaw[m]}`
                                                expr1Stat = `expr1 = ${preN} / ${preD}`;
                                                this.pyodide.runPython(`
                                                    ${expr1Stat}
                                                    if sympy.simplify(expr1) == sympy.simplify(expr):
                                                        js.flag = 1
                                                `);
                                                if (typeof(flag) !== undefined) {
                                                    return true;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            if (typeof(flag) !== undefined) {
                                return true;
                            }
                        }
                        if (typeof(flag) !== undefined) {
                            return true;
                        }
                    }
                    if (typeof(flag) !== undefined) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    /**
     * Tests whether kinetics is in the MM functional with a reactant and 2nd species as a product in the numerator
     * @param {string} kinetics string-kinetics
     * @param {[]} idsList list-id list including all the ids in kinetics, reactants and products
     * @param {[]} parametersInKineticLaw list-parameters in the kinetics
     * @param {[]} reactantList list-reactants of the reaction
     * @return {boolean}
     */
    _MMTwoSpecInNumerator(kinetics, idsList, parametersInKineticLaw, speciesInKineticLaw, reactantList) {
        var preSymbols = "";
        var i;
        for (i = 0; i < idsList.length; i++) {
            preSymbols += idsList[i];
            preSymbols += " ";
        }
        preSymbols = preSymbols.substring(0, preSymbols.length - 1);
        var preSymbolsComma = preSymbols.replaceAll(" ", ",");
        this.pyodide.runPython(`
                try:
                    ${preSymbolsComma} = sympy.symbols("${preSymbols}")
                    js.strangeFunc1 = 0
                except:
                    js.strangeFunc1 = 1
            `);
        this.pyodide.runPython(`
            try:
                expr = ${kinetics}
                js.strangeFunc2 = 0
            except:
                js.strangeFunc2 = 1
        `);
        if (strangeFunc1 + strangeFunc2 == 0) {
            if (parametersInKineticLaw.length != 0) {
                for (var j = 0; j < parametersInKineticLaw.length; j++) {
                    for (var k = 0; k < parametersInKineticLaw.length; k++) {
                        for (var l = 0; l < parametersInKineticLaw.length; l++) {
                            var cat = speciesInKineticLaw.filter(x => !reactantList.includes(x))[0];
                            var preD = `( ${reactantList[0]} + ${parametersInKineticLaw[k]} )`
                            var expr1Stat = `expr1 = ${reactantList[0]} * ${cat}  / ${preD}`
                            this.pyodide.runPython(`
                                ${expr1Stat}
                                if sympy.simplify(expr1) == sympy.simplify(expr):
                                    js.flag = 1
                            `);
                            if (typeof(flag) !== undefined) {
                                return true;
                            }
                            if (parametersInKineticLaw.length >= 2) {
                                if (k != j) {
                                    expr1Stat = `expr1 = ${reactantList[0]} * ${parametersInKineticLaw[j]} / ${preD}`;
                                    this.pyodide.runPython(`
                                        ${expr1Stat}
                                        if sympy.simplify(expr1) == sympy.simplify(expr):
                                            js.flag = 1
                                    `);
                                    if (typeof(flag) !== undefined) {
                                        return true;
                                    }
                                }
                                if (parametersInKineticLaw.length >= 3) {
                                    if (l != j && l != k) {
                                        expr1Stat = `expr1 = ${reactantList[0]} * ${parametersInKineticLaw[l]} / ${preD}`;
                                        this.pyodide.runPython(`
                                            ${expr1Stat}
                                            if sympy.simplify(expr1) == sympy.simplify(expr):
                                                js.flag = 1
                                        `);
                                        if (typeof(flag) !== undefined) {
                                            return true;
                                        }
                                    }
                                }
                            }
                        }
                        if (typeof(flag) !== undefined) {
                            return true;
                        }
                    }
                    if (typeof(flag) !== undefined) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    /**
     * Tests whether the kinetics is in the format of Hill equations.
     * 1) the numerator is one product of terms with the species to a power;
     * 2) the denomimator is the sum of two product of terms, one of which does not include the species
     *    and the other one of which include the species to the same power as the numerator.
     * @param {string} kineticsSim string-simplified kinetics
     * @param {[]} idsList list-id list including all the ids in kinetics, reactants and products
     * @param {[]} speciesInKineticLaw list-species in the kinetics
     * @return {boolean}
     */
    _HillFormat(kineticsSim, idsList, speciesInKineticLaw) {
        var flagNumerator = false;
        var flagDenominator = false;
        var eq = this._numeratorDenominator(kineticsSim, idsList);
        var numerator = eq[0];
        var denomimator = eq[1];
        var species = speciesInKineticLaw[0];

        if (!numerator.includes("+") && !numerator.includes("-")) {
            if (numerator.includes(species)) {
                if ((numerator.includes("pow(") && !numerator.includes("-1)")) || numerator.includes("**")) {
                    flagNumerator = true;
                }
            }
        }
        if (denomimator.includes("+")) {
            var terms = denomimator.split("+");
            var terms1 = terms[0];
            var terms2 = terms[1];
            if (terms1.includes(species) && !terms2.includes(species)) {
                if ((terms1.includes("pow(") && !terms1.includes("-1)")) || terms1.includes("**")) {
                    flagDenominator = true;
                }
            }
            if (terms2.includes(species) && !terms1.includes(species)) {
                if ((terms2.includes("pow(") && !terms2.includes("-1)")) || terms2.includes("**")) {
                    flagDenominator = true;
                }
            }
        }
        return flagNumerator && flagDenominator;
    }

    /**
     * Get the numerator and denominator of a "fraction" function.
     * @param {string} kineticsSim string-simplified kinetics
     * @param {[]} idsList list-id list including all the ids in kinetics, reactants and products
     * @return {[]} the numerator and the denominator of the fraction
     */
    _numeratorDenominator(kineticsSim, idsList) {
        var preSymbols = "";
        var i;
        for (i = 0; i < idsList.length; i++) {
            preSymbols += idsList[i];
            preSymbols += " ";
        }
        preSymbols = preSymbols.substring(0, preSymbols.length - 1);
        var preSymbolsComma = preSymbols.replaceAll(" ", ",");
        this.pyodide.runPython(`
                try:
                    ${preSymbolsComma} = sympy.symbols("${preSymbols}")
                    js.strangeFunc1 = 0
                except:
                    js.strangeFunc1 = 1
            `);
        this.pyodide.runPython(`
            try:
                kinetics_eq = ${kineticsSim}
                js.strangeFunc2 = 0
            except:
                js.strangeFunc2 = 1
        `);
        if (strangeFunc1 + strangeFunc2 == 0) {
            this.pyodide.runPython(`
            try: 
                js.eq0 = str(kinetics_eq.as_numer_denom()[0])
                js.eq1 = str(kinetics_eq.as_numer_denom()[1])
            except:
                js.eq0 = ""
                js.eq1 = ""
            `);
        }
        return [eq0, eq1];
    }

    /**
     * Check if a function is polynomial.
     * @param {string} kineticsSim string-simplified kinetics
     * @param {[]} idsList list-id list including all the ids in kinetics, reactants and products
     * @return {boolean}
     */
    _isPolynomial(kineticsSim, idsList) {
        if (typeof(kineticsSim))
        var preSymbols = '';
        var i;
        if (idsList.length > 0) {
            for (i = 0; i < idsList.length; i++) {
                preSymbols += idsList[i];
                preSymbols += " ";
            }
            preSymbols = preSymbols.substring(0, preSymbols.length - 1);
            var preSymbolsComma = preSymbols.replaceAll(" ", ",");
            this.pyodide.runPython(`
                try:
                    ${preSymbolsComma} = sympy.symbols("${preSymbols}")
                    js.strangeFunc1 = 0
                except:
                    js.strangeFunc1 = 1
            `);
        }

        this.pyodide.runPython(`
            try:
                kinetics_eq = ${kineticsSim}
                js.strangeFunc2 = 0
            except:
                js.strangeFunc2 = 1
        `);

        if (strangeFunc1 + strangeFunc2 == 0) {
            this.pyodide.runPython(`
                try:
                    js.polynomialFlag = kinetics_eq.is_polynomial()
                except:
                    js.polynomialFlag = false
            `);
        }
        return polynomialFlag;
    }

    _expandFormula(expansion, functionDefinitions, numRecursion, maxRecursion) {
        if (numRecursion > maxRecursion) {
            return expansion;
        }
        var done = true;
        for (var i = 0; i < functionDefinitions.length; i++) {
            var fd = functionDefinitions[i];
            var rx = new RegExp(String.raw`${fd.id}\(.*?\)`);
            var matches = rx.exec(expansion);
            var body = String(fd.body);
            if (matches == null) {
                continue;
            }
            done = false;
            for (var j = 0; j < matches.length; j++) {
                var rx_arguments = new RegExp(String.raw`\(.*?\)`);
                var argument_match = rx_arguments.exec(matches[j])[0];
                argument_match = argument_match.trim();
                argument_match = argument_match.substring(1, argument_match.length - 1);
                var call_arguments = argument_match.split(",");
                for (var k = 0; k < call_arguments.length; k++) {
                    call_arguments[k] = call_arguments[k].trim();
                }
                for (k = 0; k < fd.argumentNames.length; k++) {
                    body = body.replaceAll(fd.argumentNames[k], call_arguments[k]);
                }
                expansion = expansion.replaceAll(matches[j], body);
            }
        }
        if (!done) {
            return this._expandFormula(expansion, functionDefinitions, numRecursion + 1, maxRecursion);
        }
        return expansion;
    }

    
    _augment(astNode, result, maxDepth) {
        this._curDepth++;
        if (this._curDepth > maxDepth) {
            throw new Error("Bad Kinetics Math");
        }
        for (var i = 0; i < astNode.getNumChildren(); i++) {
            var childNode = astNode.getChild(i);
            var additions;
            if (childNode.getName() == null || childNode.getName() === '') {
                additions = this._augment(childNode, result, maxDepth);
                result = result.concat(additions);
            } else {
                if (childNode.isFunction()) {
                    additions = this._augment(childNode, result, maxDepth);
                    result = result.concat(additions);
                } else {
                    result.push(childNode.getName());
                }
            }
        }
        return result;
    }

    _getSymbols() {
        var maxDepth = 40;
        this._curDepth = 0;
        var astNode = this.libsbmlKinetics.getMath();
        var result;
        if (astNode.getName() == null || astNode.getName() === "") {
            result = [];
        } else {
            result = [astNode.getName()];
        }
        return this._augment(astNode, result, maxDepth);
    }
}

module.exports = KineticLaw;