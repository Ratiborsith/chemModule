import sqlite3

def delete_formulas_for_patent(patent_number):
    try:
        conn = sqlite3.connect('patents.db')
        cursor = conn.cursor()

        # Получаем ID патента
        cursor.execute("SELECT id FROM patents WHERE PatentNumber = ?", (patent_number,))
        patent_id = cursor.fetchone()

        if patent_id:
            patent_id = patent_id[0]

            # Удаляем все связанные записи в таблице compoundsInPatent
            cursor.execute("DELETE FROM compoundsInPatent WHERE patent_id = ?", (patent_id,))

            # Удаляем все связанные записи в таблице compounds
            cursor.execute("DELETE FROM compounds WHERE id NOT IN (SELECT compound_id FROM compoundsInPatent)")

            conn.commit()
            print(f"All formulas for patent {patent_number} have been deleted.")
        else:
            print(f"No patent found with the patent number {patent_number}.")

    except sqlite3.Error as e:
        print("An error occurred:", e)

    finally:
        if conn:
            conn.close()

# Вызываем функцию с указанным номером патента
delete_formulas_for_patent("US7754717B2")
