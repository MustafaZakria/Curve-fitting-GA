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
            totalChromosomes.add(new Chromosome(chromosome, 0));
        }
        return totalChromosomes;
    }

    private double calculateFitness(ArrayList<Double> chromosome) {

    }

    private ArrayList<Integer> calculateCumulative(ArrayList<Chromosome> chromosomes) {

        ArrayList<Integer> cumulativeFitness = new ArrayList<>();

        cumulativeFitness.add(chromosomes.get(0).getFitnessValue());
        for (int i = 1; i < chromosomes.size(); i++)
            cumulativeFitness.add(cumulativeFitness.get(i - 1) + chromosomes.get(i).getFitnessValue());

        return cumulativeFitness;

    }

    //selection and crossover repeated by 1/2 pop size
    private ArrayList<Chromosome> Selection(ArrayList<Chromosome> chromosomes) {

        ArrayList<Chromosome> selectedChromosomes = new ArrayList<>();
        ArrayList<Integer> cumulativeFitness = calculateCumulative(chromosomes);

        double randomNumber;
        int firstCumulate = cumulativeFitness.get(0);
        int lastCumulate = cumulativeFitness.get(cumulativeFitness.size() - 1);

        while (selectedChromosomes.size() != 2) {

            randomNumber = rand.nextInt(lastCumulate - firstCumulate) + firstCumulate;
            for (int i = 0; i < cumulativeFitness.size() - 1; i++) {
                if (cumulativeFitness.get(i) >= randomNumber) {
                    selectedChromosomes.add(chromosomes.get(i));
                    break;
                }
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

            int Xc = rand.nextInt(numberOfPoints - 1) + 1;

            for (int i = Xc; i < items.size(); i++) {
                int temp = afterCrossover.get(0).geneAt(i);
                afterCrossover.get(0).setGeneAt(i, afterCrossover.get(1).geneAt(i));
                afterCrossover.get(1).setGeneAt(i, temp);
            }
        }
        return afterCrossover;
    }

    private void Mutation(ArrayList<Chromosome> parents, ArrayList<Chromosome> totalChromosomes, double Pm) {

        double Rm;

        for (int i = 0; i < totalChromosomes.size(); i++) {

            for (int j = 0; j < totalChromosomes.get(i).getGenes().size(); j++) {

                Rm = rand.nextDouble() / 10;

                if (Rm <= Pm) {
                    if (totalChromosomes.get(i).geneAt(j) == 1)
                        totalChromosomes.get(i).setGeneAt(j, 0);
                    else
                        totalChromosomes.get(i).setGeneAt(j, 1);
                }
            }

            if (calculateWeight(totalChromosomes.get(i).getGenes()) <= knapsackWeight)
                totalChromosomes.get(i).setFitnessValue(calculateFitness(totalChromosomes.get(i).getGenes()));
            else {
                totalChromosomes.set(i, parents.get(i));
            }

        }
    }

    private void Replacement(ArrayList<Chromosome> oldGeneration, ArrayList<Chromosome> newGeneration) {
        oldGeneration = newGeneration;
    }

    public void solve() {

        int populationSize = 20;
        int maxGeneration = 10;
        double Pc = 0.6;
        double Pm = 0.05;
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

            Mutation(currentGeneration, newGeneration, Pm);
            Replacement(currentGeneration, newGeneration);

        }

        double bestFitness = Integer.MIN_VALUE;
        for (Chromosome c : currentGeneration) {
            if (bestFitness < c.getFitnessValue()) {
                bestFitness = c.getFitnessValue();
            }
        }
        System.out.println(bestFitness);

//        clearItems();

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
