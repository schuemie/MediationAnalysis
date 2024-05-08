CREATE TABLE #Codesets (
  codeset_id int NOT NULL,
  concept_id bigint NOT NULL
)
;

INSERT INTO #Codesets (codeset_id, concept_id)
SELECT 0 as codeset_id, c.concept_id FROM (select distinct I.concept_id FROM
( 
  select concept_id from @vocabulary_database_schema.CONCEPT where concept_id in (4291005,4029488,4112183,24966,46273476,4001172,4340383,4055225,4095432,4003021,4340948,46269836,4012113,4159144,28779,4026125,201519,4059284,4058695,194417,197917,4240725,192670,4059297,201613,439673,40490929,4109883,4308946,4098652,200762,197795,4340385,192675,4059299,22340,4313846,193256,4001664,4310665,192680,45769824,4245975,194856,4058694,197493,4340941,4111998,40484532,4046123,198683,192242,196455,439675,439674,4059290,196463,201612,46269835,4026136,193693,4246127,4064161,4238978,4001171,46269816,200031,4340386,45768827,201343,200763,4135822,46269818,192240,40481531,4277276,201901,763021,194087,4281232,4340394,46270037,4055224,439672,196625,194984,196029,200451,194692,199867,197494,4059298,194990,4340390,4115271,4092685,4267417,4058696,198964,197490,4211974)

) I
) C UNION ALL 
SELECT 1 as codeset_id, c.concept_id FROM (select distinct I.concept_id FROM
( 
  select concept_id from @vocabulary_database_schema.CONCEPT where concept_id in (2735189,2735471,2735462,2735468,2735460,2735466,2735461,2735467,2735464,2735470,2735463,2735469,2735214,2735220,2735211,2735217,2735209,2735215,2735210,2735216,2735213,2735219,2735212,2735218,2735429,2735435,2735426,2735432,2735424,2735430,2735425,2735431,2735428,2735434,2735427,2735433,2735453,2735459,2735450,2735456,2735448,2735454,2735449,2735455,2735452,2735458,2735451,2735457,2734964,2735171,2734955,2734973,2734949,2734967,2734952,2734970,2734961,2735168,2734958,2734976,2734965,2735172,2734956,2734974,2734950,2734968,2734953,2734971,2734962,2735169,2734959,2735166,2734963,2735170,2734954,2734972,2734948,2734966,2734951,2734969,2734960,2735167,2734957,2734975,2735980,2735986,2735977,2735983,2735975,2735981,2735976,2735982,2735979,2735985,2735978,2735984,2735718,2735724,2735715,2735721,2735713,2735719,2735714,2735720,2735717,2735723,2735716,2735722,2735679,2735699,2735479,2735689,2735473,2735683,2735476,2735686,2735676,2735696,2735673,2735693,2735680,2735700,2735480,2735690,2735474,2735684,2735477,2735687,2735681,2735691,2735677,2735697,2735674,2735694,2735678,2735698,2735478,2735688,2735472,2735682,2735475,2735685,2735675,2735695,2735481,2735692,2735706,2735712,2735703,2735709,2735701,2735707,2735702,2735708,2735705,2735711,2735704,2735710,4080760,4284254,2735465,2735207,2735180,2735198,2735174,2735192,2735177,2735195,2735186,2735204,2735183,2735201,2735190,2735208,2735181,2735199,2735175,2735193,2735178,2735196,2735187,2735205,2735184,2735202,2735188,2735206,2735179,2735197,2735173,2735191,2735176,2735194,2735185,2735203,2735182,2735200,2735441,2735447,2735438,2735444,2735436,2735442,2735437,2735443,2735440,2735446,2735439,2735445,2740046,2740049,2740052)

) I
) C
;

with primary_events (event_id, person_id, start_date, end_date, op_start_date, op_end_date, visit_occurrence_id) as
(
-- Begin Primary Events
select P.ordinal as event_id, P.person_id, P.start_date, P.end_date, op_start_date, op_end_date, cast(P.visit_occurrence_id as bigint) as visit_occurrence_id
FROM
(
  select E.person_id, E.start_date, E.end_date,
         row_number() OVER (PARTITION BY E.person_id ORDER BY E.sort_date ASC) ordinal,
         OP.observation_period_start_date as op_start_date, OP.observation_period_end_date as op_end_date, cast(E.visit_occurrence_id as bigint) as visit_occurrence_id
  FROM 
  (
  -- Begin Condition Occurrence Criteria
SELECT C.person_id, C.condition_occurrence_id as event_id, C.condition_start_date as start_date, COALESCE(C.condition_end_date, DATEADD(day,1,C.condition_start_date)) as end_date,
  C.visit_occurrence_id, C.condition_start_date as sort_date
FROM 
(
  SELECT co.* 
  FROM @cdm_database_schema.CONDITION_OCCURRENCE co
  JOIN #Codesets codesets on ((co.condition_concept_id = codesets.concept_id and codesets.codeset_id = 0))
) C


-- End Condition Occurrence Criteria

UNION ALL
-- Begin Procedure Occurrence Criteria
select C.person_id, C.procedure_occurrence_id as event_id, C.procedure_date as start_date, DATEADD(d,1,C.procedure_date) as END_DATE,
       C.visit_occurrence_id, C.procedure_date as sort_date
from 
(
  select po.* 
  FROM @cdm_database_schema.PROCEDURE_OCCURRENCE po
JOIN #Codesets codesets on ((po.procedure_concept_id = codesets.concept_id and codesets.codeset_id = 1))
) C


-- End Procedure Occurrence Criteria

  ) E
	JOIN @cdm_database_schema.observation_period OP on E.person_id = OP.person_id and E.start_date >=  OP.observation_period_start_date and E.start_date <= op.observation_period_end_date
  WHERE DATEADD(day,0,OP.OBSERVATION_PERIOD_START_DATE) <= E.START_DATE AND DATEADD(day,0,E.START_DATE) <= OP.OBSERVATION_PERIOD_END_DATE
) P
WHERE P.ordinal = 1
-- End Primary Events

)
SELECT event_id, person_id, start_date, end_date, op_start_date, op_end_date, visit_occurrence_id
INTO #qualified_events
FROM 
(
  select pe.event_id, pe.person_id, pe.start_date, pe.end_date, pe.op_start_date, pe.op_end_date, row_number() over (partition by pe.person_id order by pe.start_date ASC) as ordinal, cast(pe.visit_occurrence_id as bigint) as visit_occurrence_id
  FROM primary_events pe
  
) QE

;

--- Inclusion Rule Inserts

create table #inclusion_events (inclusion_rule_id bigint,
	person_id bigint,
	event_id bigint
);

with cteIncludedEvents(event_id, person_id, start_date, end_date, op_start_date, op_end_date, ordinal) as
(
  SELECT event_id, person_id, start_date, end_date, op_start_date, op_end_date, row_number() over (partition by person_id order by start_date ASC) as ordinal
  from
  (
    select Q.event_id, Q.person_id, Q.start_date, Q.end_date, Q.op_start_date, Q.op_end_date, SUM(coalesce(POWER(cast(2 as bigint), I.inclusion_rule_id), 0)) as inclusion_rule_mask
    from #qualified_events Q
    LEFT JOIN #inclusion_events I on I.person_id = Q.person_id and I.event_id = Q.event_id
    GROUP BY Q.event_id, Q.person_id, Q.start_date, Q.end_date, Q.op_start_date, Q.op_end_date
  ) MG -- matching groups

)
select event_id, person_id, start_date, end_date, op_start_date, op_end_date
into #included_events
FROM cteIncludedEvents Results
WHERE Results.ordinal = 1
;



-- generate cohort periods into #final_cohort
with cohort_ends (event_id, person_id, end_date) as
(
	-- cohort exit dates
  -- By default, cohort exit at the event's op end date
select event_id, person_id, op_end_date as end_date from #included_events
),
first_ends (person_id, start_date, end_date) as
(
	select F.person_id, F.start_date, F.end_date
	FROM (
	  select I.event_id, I.person_id, I.start_date, E.end_date, row_number() over (partition by I.person_id, I.event_id order by E.end_date) as ordinal 
	  from #included_events I
	  join cohort_ends E on I.event_id = E.event_id and I.person_id = E.person_id and E.end_date >= I.start_date
	) F
	WHERE F.ordinal = 1
)
select person_id, start_date, end_date
INTO #cohort_rows
from first_ends;

with cteEndDates (person_id, end_date) AS -- the magic
(	
	SELECT
		person_id
		, DATEADD(day,-1 * 0, event_date)  as end_date
	FROM
	(
		SELECT
			person_id
			, event_date
			, event_type
			, MAX(start_ordinal) OVER (PARTITION BY person_id ORDER BY event_date, event_type ROWS UNBOUNDED PRECEDING) AS start_ordinal 
			, ROW_NUMBER() OVER (PARTITION BY person_id ORDER BY event_date, event_type) AS overall_ord
		FROM
		(
			SELECT
				person_id
				, start_date AS event_date
				, -1 AS event_type
				, ROW_NUMBER() OVER (PARTITION BY person_id ORDER BY start_date) AS start_ordinal
			FROM #cohort_rows
		
			UNION ALL
		

			SELECT
				person_id
				, DATEADD(day,0,end_date) as end_date
				, 1 AS event_type
				, NULL
			FROM #cohort_rows
		) RAWDATA
	) e
	WHERE (2 * e.start_ordinal) - e.overall_ord = 0
),
cteEnds (person_id, start_date, end_date) AS
(
	SELECT
		 c.person_id
		, c.start_date
		, MIN(e.end_date) AS end_date
	FROM #cohort_rows c
	JOIN cteEndDates e ON c.person_id = e.person_id AND e.end_date >= c.start_date
	GROUP BY c.person_id, c.start_date
)
select person_id, min(start_date) as start_date, end_date
into #final_cohort
from cteEnds
group by person_id, end_date
;

DELETE FROM @target_database_schema.@target_cohort_table where cohort_definition_id = @target_cohort_id;
INSERT INTO @target_database_schema.@target_cohort_table (cohort_definition_id, subject_id, cohort_start_date, cohort_end_date)
select @target_cohort_id as cohort_definition_id, person_id, start_date, end_date 
FROM #final_cohort CO
;





TRUNCATE TABLE #cohort_rows;
DROP TABLE #cohort_rows;

TRUNCATE TABLE #final_cohort;
DROP TABLE #final_cohort;

TRUNCATE TABLE #inclusion_events;
DROP TABLE #inclusion_events;

TRUNCATE TABLE #qualified_events;
DROP TABLE #qualified_events;

TRUNCATE TABLE #included_events;
DROP TABLE #included_events;

TRUNCATE TABLE #Codesets;
DROP TABLE #Codesets;

