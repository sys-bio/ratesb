/** Provides Information on SBML Kinetics Laws */

class KineticLaw {
    constructor( libsbmlModel, libsbmlReaction, libsbmlKinetics, reaction, functionDefinitions, pyodide, processSBML) {
        this.libsbmlModel = libsbmlModel;
        this.libsbmlReaction = libsbmlReaction;
        this.pyodide = pyodide;
        this.libsbmlKinetics = libsbmlKinetics;
        console.log(this.libsbmlKinetics);
        this.formula = this.libsbmlKinetics.getFormula();
        this.reaction = reaction;
        this.model = processSBML;
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
        this.analysis = {};
        this.classify();
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
        var parametersInKineticLaw = [];
        var othersInKineticLaw = [];
        var idsList = [...new Set(this.symbols)];
        var i;
        for (i = 0; i < idsList.length; i++) {
            if (this.model.speciesList.includes(idsList[i])) {
                speciesInKineticLaw.push(idsList[i]);
            } else if (this.model.parameterList.includes(idsList[i])) {
                parametersInKineticLaw.push(idsList[i]);
            } else {
                othersInKineticLaw.push(idsList[i]);
            }
        }
        parametersInKineticLaw = parametersInKineticLaw.concat(othersInKineticLaw);
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
        this.classificationCp["MMcat"] = this.isMMcat(kinetics, kineticsSim, idsList, speciesInKineticLaw, parametersInKineticLaw, reactantList);
        this.classificationCp["Hill"] = this.isHill(kineticsSim, idsList, speciesInKineticLaw);
        this.classificationCp["Fraction"] = this.isFraction(kineticsSim, idsList, speciesInKineticLaw);
        this.classificationCp["Polynomial"] = this.isPolynomial(kineticsSim, idsList, speciesInKineticLaw);

        console.log(this.classificationCp);
        this.kineticLawAnalysis(speciesInKineticLaw, parametersInKineticLaw, othersInKineticLaw, idsList, kinetics, kineticsSim, reactantList, productList);
    }

    kineticLawAnalysis(speciesInKineticLaw, parametersInKineticLaw, othersInKineticLaw, idsList, kinetics, kineticsSim, reactantList, productList) {
        this.analysis["floatingSpecies"] = this.checkFloatingSpecies(speciesInKineticLaw, reactantList, productList);
        this.analysis["nonIncreasingSpecies"] = this.checkNonIncreasingSpecies(speciesInKineticLaw, parametersInKineticLaw, reactantList, productList, kinetics);
        console.log(this.analysis);
    }

    /**
     * Check if the rate law contains all reactants
     * If reversible, check if the rate law contains all products
     * Return a list of species that violates this rule
     * @param {[]} speciesInKineticLaw list-species in the kinetics
     * @param {[]} reactantList list-reactants of the reaction
     * @param {[]} productList list-products of the reaction
     * @return {[]} list of species violating the rule
     */
    checkFloatingSpecies(speciesInKineticLaw, reactantList, productList) {
        var floatingSpecies = [];
        reactantList.forEach(reactant => {
            if (!speciesInKineticLaw.includes(reactant)) {
                floatingSpecies.push(reactant);
            }
        });
        if (this.libsbmlReaction.getReversible()) {
            productList.forEach(product => {
                if (!speciesInKineticLaw.includes(product)) {
                    floatingSpecies.push(product);
                }
            });
        }
        return floatingSpecies;
    }

    /**
     * Check if the rate law is increasing as each reactant increases
     * If reversible, check if the rate law is decreasing as each product increases
     * Return a list of species that violates this rule
     * @param {[]} speciesInKineticLaw list-species in the kinetics
     * @param {[]} parametersInKineticLaw list-parameters in the kinetics
     * @param {[]} reactantList list-reactants of the reaction
     * @param {[]} productList list-products of the reaction
     * @param {string} kinetics string-kinetics
     * @return {[]} list of species violating the rule
     */
    checkNonIncreasingSpecies(speciesInKineticLaw, parametersInKineticLaw, reactantList, productList, kinetics) {
        console.log(speciesInKineticLaw)
        console.log(parametersInKineticLaw)
        console.log(reactantList)
        console.log(productList)
        console.log(kinetics)
        var nonIncreasingSpecies = [];
        var i;
        for (i = 0; i < this.libsbmlModel.getNumCompartments(); i++) {
            if (parametersInKineticLaw.includes(this.libsbmlModel.getCompartment(i).getId())) {
                console.log(this.libsbmlModel.getCompartment(i))
                this.pyodide.runPython(`
                    ${this.libsbmlModel.getCompartment(i).getId()} = 1
                `);
            }
        }
        var param;
        for (i = 0; i < this.libsbmlKinetics.getNumParameters(); i++) {
            param = this.libsbmlKinetics.getParameter(i);
            this.pyodide.runPython(`
                ${param.getId()} =  ${param.getValue() > 0 ? 1 : param.getValue() < 0 ? -1 : 0}
            `);
        }
        for (i = 0; i < this.libsbmlKinetics.getNumLocalParameters(); i++) {
            param = this.libsbmlKinetics.getLocalParameter(i);
            this.pyodide.runPython(`
                ${param.getId()} = ${param.getValue() > 0 ? 1 : param.getValue() < 0 ? -1 : 0}
            `);
        }
        reactantList.forEach(reactant => {
            speciesInKineticLaw.forEach(otherSpecies => {
                if (otherSpecies !== reactant) {
                    this.pyodide.runPython(`
                        ${otherSpecies} = 1
                    `);
                }
            });
            console.log(reactant)
            this.pyodide.runPython(`
                diff_kinetics = sympy.diff(${kinetics}, ${reactant})
                ddx = sympy.lambdify(${reactant}, diff_kinetics)
                js.isIncreasing = ddx(1) > 0
            `);
            if (!isIncreasing) {
                nonIncreasingSpecies.push(reactant);
            }
            speciesInKineticLaw.forEach(otherSpecies => {
                this.pyodide.runPython(`
                    ${otherSpecies} = sympy.symbols("${otherSpecies}")
                `);
            });
        });
        if (this.libsbmlReaction.getReversible()) {
            productList.forEach(product => {
                productList.forEach(otherProduct => {
                    if (otherProduct !== product) {
                        this.pyodide.runPython(`
                            ${otherProduct} = 1
                        `);
                    }
                });
                this.pyodide.runPython(`
                    diff_kinetics = sympy.diff(${kinetics}, ${product})
                    ddu = sympy.lambdify(${product}, diff_kinetics)
                    js.isIncreasing = ddx(1) < 0
                `);
                if (!isIncreasing) {
                    nonIncreasingSpecies.push(product);
                }
                speciesInKineticLaw.forEach(otherSpecies => {
                    this.pyodide.runPython(`
                        ${otherSpecies} = sympy.symbols("${otherSpecies}")
                    `);
                });
            });
        }
        return nonIncreasingSpecies;
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

    isMMcat(kinetics, kineticsSim, idsList, speciesInKineticLaw, parametersInKineticLaw, reactantList) {
        var eq = this._numeratorDenominator(kineticsSim, idsList);
        var flagFr = false;
        speciesInKineticLaw.forEach(species => {
            if (eq[1].includes(species)) {
                flagFr = true;
            }
        });
        if (flagFr) {
            if (this._numSpeciesInKinetics(speciesInKineticLaw) == 2 && this._numOfRcts(reactantList) == 1) {
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
                if ((!numerator.includes("pow(") && !numerator.includes("-1)")) || numerator.includes("**")) {
                    flagNumerator = true;
                }
            }
        }
        if (denomimator.includes("+")) {
            var terms = denomimator.split("+");
            var terms1 = terms[0];
            var terms2 = terms[1];
            if (terms1.includes(species) && !terms2.includes(species)) {
                if ((!terms1.includes("pow(") && !terms1.includes("-1)")) || terms1.includes("**")) {
                    flagDenominator = true;
                }
            }
            if (terms2.includes(species) && !terms1.includes(species)) {
                if ((!terms2.includes("pow(") && !terms2.includes("-1)")) || terms2.includes("**")) {
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