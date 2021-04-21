package online.operators;

import beast.core.Description;
import beast.core.Operator;
import beast.core.OperatorSchedule;

@Description("Operator schedule that initially only uses proposals from PartitionOperators "
		+ "then later include all operators")
public class AfterburnOperatorSchedule extends OperatorSchedule {

	
	private long count = 0;
	private long limit = 1000;
	
	public void reset(long limit) {
		count = 0;
		this.limit = limit;
	}
	
	@Override
	public Operator selectOperator() {
		count++;
		if (count < limit) {
			int attempt = 0;
			while (attempt < 1000) {
				Operator operator = super.selectOperator();
				if (operator instanceof PartitionOperator) {
					return operator;
				}
			}
			// there are probably no PartitionOperators
			// carry on with a standard operator
		}
		if (subschedulesInput.get().size() != 0) {
			return subschedulesInput.get().get(0).selectOperator();
		}
		return super.selectOperator();
	}
}
