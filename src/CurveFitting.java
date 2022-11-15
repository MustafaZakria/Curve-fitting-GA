import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;

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

    public static void main(String[] args) {

        Scanner sc = new Scanner(System.in);
        int testCases = sc.nextInt();

        int numberOfPoints, degree, x, y;
        ArrayList<Point> points = new ArrayList<>();

        for (int i = 0; i < testCases; i++) {

            numberOfPoints = sc.nextInt();
            degree = sc.nextInt();

            for (int j = 0; j < numberOfPoints; j++) {
                x = sc.nextInt();
                y = sc.nextInt();
                points.add(new Point(x, y));
            }

            CurveFitting cf = new CurveFitting(degree, numberOfPoints, points);
            cf.solve();

        }

    }

    private double calculateError(ArrayList<Double> chromosome) {

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

        return 1.0 / error;

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

            totalChromosomes.add(new Chromosome(chromosome, calculateError(chromosome)));

        }

        return totalChromosomes;

    }

    //selection and crossover repeated by 1/2 pop size
    private ArrayList<Chromosome> Selection(ArrayList<Chromosome> chromosomes) {

        ArrayList<Chromosome> selectedChromosomes = new ArrayList<>();

        int randomNumber;

        while (selectedChromosomes.size() != 2) {

            randomNumber = rand.nextInt(chromosomes.size() + 1);
            selectedChromosomes.add(chromosomes.get(randomNumber));

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

    private void Mutation(ArrayList<Chromosome> totalChromosomes, double Pm, int currentGeneration, int maxGeneration, int dependencyFactor) {

        double Rm;
        double lowerBound, upperBound, deltaLower, deltaUpper, y, delta;
        double random;

        for (Chromosome chromosome : totalChromosomes) {

            lowerBound = getLowerBound(chromosome);
            upperBound = getUpperBound(chromosome);

            Rm = rand.nextDouble() / 10;

            for (int i = 0; i < chromosome.getGenes().size(); i++) {

                if (Rm <= Pm) {

                    random = rand.nextDouble();
                    deltaLower = chromosome.getGeneAt(i) - lowerBound;
                    deltaUpper = upperBound - chromosome.getGeneAt(i);

                    if (random <= 0.5)
                        y = deltaLower;
                    else
                        y = deltaUpper;

                    delta = y * (Math.pow(1 - random, Math.pow(1 - (double)currentGeneration/maxGeneration, dependencyFactor)));

                    if (y == deltaLower)
                        chromosome.setGeneAt(i, chromosome.getGeneAt(i) - delta);
                    else if (y == deltaUpper)
                        chromosome.setGeneAt(i, chromosome.getGeneAt(i) + delta);


                }

            }

        }

    }

    private void Replacement(ArrayList<Chromosome> oldGeneration, ArrayList<Chromosome> newGeneration) {

    }

    public void solve() {

        int populationSize = 20;
        int maxGeneration = 10;
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
            Replacement(currentGeneration, newGeneration);

        }

//        double bestFitness = Integer.MIN_VALUE;
//        for (Chromosome c : currentGeneration) {
//            if (bestFitness < c.getFitnessValue()) {
//                bestFitness = c.getFitnessValue();
//            }
//        }
//        System.out.println(bestFitness);

        clearPoints();

    }

    static class Chromosome {

        private ArrayList<Double> genes;
        private double fitnessValue;

        public Chromosome(ArrayList<Double> genes, double fitnessValue) {
            this.genes = genes;
            this.fitnessValue = fitnessValue;
        }

        public ArrayList<Double> getGenes() {
            return genes;
        }

        public void setGenes(ArrayList<Double> genes) {
            this.genes = genes;
        }

        public double getFitnessValue() {
            return fitnessValue;
        }

        public void setFitnessValue(double fitnessValue) {
            this.fitnessValue = fitnessValue;
        }

        public double getGeneAt(int pos) {
            return genes.get(pos);
        }

        public void setGeneAt(int pos, double value) {
            genes.set(pos, value);
        }

    }

    static class Point {
        private int x;
        private int y;

        public Point(int x, int y) {
            this.x = x;
            this.y = y;
        }

        public void setX(int x) {
            this.x = x;
        }

        public void setY(int y) {
            this.y = y;
        }

        public int getX() {
            return x;
        }

        public int getY() {
            return y;
        }


    }

}
