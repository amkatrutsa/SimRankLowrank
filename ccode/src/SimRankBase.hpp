
class SimRankBase {
    public:
        SimRankBase() {}
        explicit SimRankBase(double c, int num_iter) : c_(c), num_iter_(num_iter) {}
        virtual bool compute() = 0;
    protected:
        double c_;
        int num_iter_;
};