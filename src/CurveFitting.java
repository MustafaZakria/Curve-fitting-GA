import java.awt.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.util.*;

public class CurveFitting {

    Random rand = new Random();
    private int numberOfPoints;
    private int degree;
    private ArrayList<Point> points;

    public CurveFitting(int degree, int numberOfPoints, ArrayList<Point> points) {
        this.degree = degree;
        this.numberOfPoints = numberOfPoints;
        this.points = points;
    }

    public int getNumberOfPoints() {
        return numberOfPoints;
    }

    public void setNumberOfPoints(int numberOfPoints) {
        this.numberOfPoints = numberOfPoints;
    }

    public int getDegree() {
        return degree;
    }

    public void setDegree(int degree) {
        this.degree = degree;
    }

    public ArrayList<Point> getPoints() {
        return points;
    }

    public void setPoints(ArrayList<Point> points) {
        this.points = points;
    }

    public void clearPoints() {
        points.clear();
    }

    public static void main(String[] args) throws FileNotFoundException {

        Scanner scan;
        String result="";
        File inputFile = new File("C:\\Users\\lenovo\\Desktop\\Soft computing\\Assign2-Curve-fitting-GA\\input.txt");
        try {
            scan = new Scanner(inputFile);

            while (scan.hasNextDouble()) {
                double testCases = scan.nextDouble();

                double numberOfPoints, degree;
                double x, y;
                ArrayList<Point> points = new ArrayList<>();

                for (int i = 0; i < testCases; i++) {

                    numberOfPoints = scan.nextDouble();
                    degree = scan.nextDouble();

                    for (int j = 0; j < numberOfPoints; j++) {
                        x = scan.nextDouble();
                        y = scan.nextDouble();
                        points.add(new Point(x, y));
                    }

                    CurveFitting cf = new CurveFitting((int) degree, (int) numberOfPoints, points);
                    result += i + "\n" + cf.solve() + "\n" ;
                }
                writeToFile(result);

            }
        } catch (FileNotFoundException e1) {
            e1.printStackTrace();
        }
    }

    private static void writeToFile(String string) {

        File outputFile = new File("output.txt");
        try {
            outputFile.createNewFile();
            FileWriter myWriter = new FileWriter("output.txt");
            myWriter.write(string);
            myWriter.close();
        }
        catch(Exception e) {
            e.getStackTrace();
        }
    }

    public String solve() {

        int populationSize = 80;
        int maxGeneration = 50;
        int dependencyFactor = 2;
        double Pc = 0.6;
        double Pm = 0.02;

        ArrayList<Chromosome> currentGeneration = InitializePopulation(populationSize);

        // INSIDE A FOR LOOP FOR MAX GENERATIONS
        for (int i = 0; i < maxGeneration; i++) {

            ArrayList<Chromosome> newGeneration = new ArrayList<>();

            for (int j = 0; j < populationSize / 2; j++) {
                ArrayList<Chromosome> selectedParents = Selection(currentGeneration);
                ArrayList<Chromosome> parentsCrossover = Crossover(selectedParents, Pc);
                newGeneration.add(parentsCrossover.get(0));
                newGeneration.add(parentsCrossover.get(1));
            }

            Mutation(newGeneration, Pm, i, maxGeneration, dependencyFactor);
            Replacement(currentGeneration, newGeneration, populationSize);

        }

        double lowestError = Double.MAX_VALUE;
        Chromosome bestChromosome = new Chromosome();
        for (Chromosome c : currentGeneration) {
            if (lowestError > c.getError()) {
                lowestError = c.getError();
                bestChromosome = c;
            }
        }

        clearPoints();
        String result = bestChromosome.getGenes().toString() + "\t" + lowestError;
        return result;

    }


    private ArrayList<Chromosome> InitializePopulation(int populationSize) {

        ArrayList<Chromosome> totalChromosomes = new ArrayList<>(populationSize);

        int chromosomeSize = degree + 1;
        int rangeMin = -10;
        int rangeMax = 10;

        for (int i = 0; i < populationSize; i++) {

            ArrayList<Double> chromosome = new ArrayList<>(chromosomeSize);
            double random;

            for (int gene = 0; gene < chromosomeSize; gene++) {
                random = rangeMin + (rangeMax - rangeMin) * rand.nextDouble();
                chromosome.add(random);
            }

            totalChromosomes.add(new Chromosome(chromosome, calculateFitness(chromosome)));

        }

        return totalChromosomes;

    }

    //selection and crossover repeated by 1/2 pop size
    private ArrayList<Chromosome> Selection(ArrayList<Chromosome> chromosomes) {

        ArrayList<Chromosome> selectedChromosomes = new ArrayList<>();
        ArrayList<Chromosome> compareTwoChromosomes = new ArrayList<>();

        int randomNumber;

        while (selectedChromosomes.size() != 2) {

            randomNumber = rand.nextInt(chromosomes.size());
            compareTwoChromosomes.add(chromosomes.get(randomNumber));
            if (compareTwoChromosomes.size() == 2) {
                if (compareTwoChromosomes.get(0).getError() < compareTwoChromosomes.get(1).getError())
                    selectedChromosomes.add(compareTwoChromosomes.get(0));
                else
                    selectedChromosomes.add(compareTwoChromosomes.get(1));
                compareTwoChromosomes.clear();
            }
        }

        return selectedChromosomes;

    }

    private ArrayList<Chromosome> Crossover(ArrayList<Chromosome> selectedChromosomes, double Pc) {

        if (selectedChromosomes.get(0) == selectedChromosomes.get(1))
            return selectedChromosomes;

        double Rc = rand.nextDouble();

        ArrayList<Chromosome> afterCrossover = new ArrayList<>();
        afterCrossover.add(selectedChromosomes.get(0));
        afterCrossover.add(selectedChromosomes.get(1));


        if (Rc <= Pc) {

            int firstPoint = rand.nextInt(degree) + 1;
            int secondPoint = firstPoint;

            while (firstPoint == secondPoint) {
                secondPoint = rand.nextInt(degree) + 1;
            }

            if (firstPoint > secondPoint) {
                int temp = firstPoint;
                firstPoint = secondPoint;
                secondPoint = temp;
            }

            Chromosome tempSelected = afterCrossover.get(0);

            for (int i = firstPoint; i < secondPoint; i++) {
                afterCrossover.get(0).setGeneAt(i, afterCrossover.get(1).getGeneAt(i));
                afterCrossover.get(1).setGeneAt(i, tempSelected.getGeneAt(i));
            }

        }

        return afterCrossover;

    }


    private void Mutation(ArrayList<Chromosome> totalChromosomes, double Pm, int currentGeneration, int maxGeneration, int dependencyFactor) {

        double Rm;
        double lowerBound, upperBound, deltaLower, deltaUpper, y, delta;
        double random;

        for (Chromosome chromosome : totalChromosomes) {

            lowerBound = getLowerBound(chromosome);
            upperBound = getUpperBound(chromosome);

            for (int i = 0; i < chromosome.getGenes().size(); i++) {

                Rm = rand.nextDouble() / 10;

                if (Rm <= Pm) {

                    deltaLower = chromosome.getGeneAt(i) - lowerBound;
                    deltaUpper = upperBound - chromosome.getGeneAt(i);
                    random = rand.nextDouble();

                    if (random <= 0.5)
                        y = deltaLower;
                    else
                        y = deltaUpper;

                    delta = y * (Math.pow(1 - random, Math.pow(1 - (double) currentGeneration / maxGeneration, dependencyFactor)));

                    if (y == deltaLower)
                        chromosome.setGeneAt(i, chromosome.getGeneAt(i) - delta);
                    else if (y == deltaUpper)
                        chromosome.setGeneAt(i, chromosome.getGeneAt(i) + delta);

                }

            }

        }

    }

    private void Replacement(ArrayList<Chromosome> oldGeneration, ArrayList<Chromosome> newGeneration, int popSize) {
        ArrayList<Chromosome> twoGenerations = new ArrayList<>();
        twoGenerations.addAll(oldGeneration);
        twoGenerations.addAll(newGeneration);
        twoGenerations.sort((o1, o2) -> Double.compare(o2.getError(), o1.getError()));
        for (int i = 0; i < popSize; i++) {
            oldGeneration.set(i, twoGenerations.get(twoGenerations.size() - 1 - i));
//            System.out.println(oldGeneration.get(i).getError());
        }
    }

    private double calculateFitness(ArrayList<Double> chromosome) {

        double error = 0.0;

        for (int i = 0; i < numberOfPoints; i++) {

            double fitness = 0.0;
            for (int j = 0; j < chromosome.size(); j++) {
                fitness += chromosome.get(j) * Math.pow(points.get(i).getX(), j);
            }
            fitness -= points.get(i).getY();
            fitness = Math.pow(fitness, 2);

            error += fitness;

        }
        error /= numberOfPoints;
        //return 1.0 / error;
        return error;  //the fittest is the lowest fitness value
    }

    private double getLowerBound(Chromosome chromosome) {

        double lower = Double.MAX_VALUE;
        for (int i = 0; i < chromosome.getGenes().size(); i++) {
            if (lower > chromosome.getGeneAt(i)) {
                lower = chromosome.getGeneAt(i);
            }
        }
        return lower;

    }

    private double getUpperBound(Chromosome chromosome) {

        double upper = Double.MIN_VALUE;
        for (int i = 0; i < chromosome.getGenes().size(); i++) {
            if (upper < chromosome.getGeneAt(i)) {
                upper = chromosome.getGeneAt(i);
            }
        }
        return upper;

    }

    static class Chromosome {

        private ArrayList<Double> genes;
        private double error;

        public Chromosome() {
            genes = new ArrayList<>();
            error = 0.0;
        }

        public Chromosome(ArrayList<Double> genes, double error) {
            this.genes = genes;
            this.error = error;
        }

        public ArrayList<Double> getGenes() {
            return genes;
        }

        public void setGenes(ArrayList<Double> genes) {
            this.genes = genes;
        }

        public double getError() {
            return error;
        }

        public void setFitnessValue(double error) {
            this.error = error;
        }

        public double getGeneAt(int pos) {
            return genes.get(pos);
        }

        public void setGeneAt(int pos, double value) {
            genes.set(pos, value);
        }

    }

    static class Point {
        private double x;
        private double y;

        public Point(double x, double y) {
            this.x = x;
            this.y = y;
        }

        public double getX() {
            return x;
        }

        public void setX(double x) {
            this.x = x;
        }

        public double getY() {
            return y;
        }

        public void setY(double y) {
            this.y = y;
        }
    }

}
